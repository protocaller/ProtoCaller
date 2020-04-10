import copy as _copy
import logging as _logging
import os as _os
import tarfile as _tarfile
import time as _time
import warnings as _warnings

from selenium.common import exceptions as _exceptions
from selenium.webdriver.firefox import options as _options
from selenium.webdriver.common import by as _by
from selenium.webdriver.support import expected_conditions as _EC
from selenium.webdriver.support import wait as _wait
import seleniumrequests as _seleniumrequests

import ProtoCaller as _PC
from ProtoCaller.IO.PDB import PDB as _PDB

__all__ = ["PDBReader", "ligandReader", "charmmguiTransform"]


def PDBReader(file, timeout=600):
    """
    Accesses http://charmm-gui.org and uses the PDB Reader.

    Parameters
    ----------
    file : str
        Path to the input PDB file.
    timeout : int
        Timeout in seconds.

    Returns
    -------
    filename_output : str
        The absolute path to the output TGZ archive.
    """
    def autoClicker(id, timeout):
        # deals with some rare cases of an unclickable element
        for i in range(timeout):
            try:
                elem = driver.find_element_by_id(id)
                elem.click()
                return
            except _exceptions.WebDriverException:
                _time.sleep(1)
        elem = driver.find_element_by_id(id)
        elem.click()

    file = _os.path.abspath(file)
    options = _options.Options()
    options.headless = _PC.HEADLESS_CHARMMGUI

    try:
        driver = _seleniumrequests.Chrome(options=options)
    except _exceptions.WebDriverException:
        try:
            driver = _seleniumrequests.Firefox(options=options)
        except _exceptions.WebDriverException:
            raise SystemError("Need either Chrome or Firefox for CHARMM-GUI "
                              "functionality.")

    _logging.info("Accessing http://www.charmm-gui.org ...")
    driver.get("http://www.charmm-gui.org/?doc=input/pdbreader")

    pdb_element = driver.find_element_by_name("file")
    pdb_element.send_keys(file)

    pdb_radio = driver.find_element_by_xpath("//input[@name='pdb_format' and "
                                             "@value='PDB']")
    pdb_radio.click()

    autoClicker("nextBtn", 60)

    # could add some support for options. For now, we just go with the
    # defaults.
    wait = _wait.WebDriverWait(driver, timeout)
    wait.until(lambda driver: driver.current_url == "http://www.charmm-gui.org/?doc=input/pdbreader&step=1")
    wait.until(_EC.element_to_be_clickable((_by.By.ID, "nextBtn")))
    autoClicker("nextBtn", 60)

    wait.until(lambda driver: driver.current_url == "http://www.charmm-gui.org/?doc=input/pdbreader&step=2")
    wait.until(_EC.element_to_be_clickable((_by.By.ID, "nextBtn")))
    autoClicker("nextBtn", 60)

    try:
        _logging.info("Retrieving files...")
        wait.until(lambda driver: driver.current_url == "http://www.charmm-gui.org/?doc=input/pdbreader&step=3")
        wait.until(_EC.visibility_of_any_elements_located((_by.By.CLASS_NAME,
                                                           "download")))
    except TimeoutError:
        raise ConnectionError("Could not retrieve any files. Please increase "
                              "the maximum timeout or try again later.")

    _logging.info("Downloading TGZ archive...")
    filebase = _os.path.splitext(file)[0]
    tgz_file = driver.find_elements_by_partial_link_text(".tgz")[0]
    response = driver.request('POST', tgz_file.get_attribute("href"),
                              verify=False, stream=True)
    with open(filebase + "_CHARMM.tgz", "wb") as file:
        file.write(response.raw.read())

    driver.quit()
    return filebase + "_CHARMM.tgz"


def ligandReader(file, timeout=60, find_similar_residues=False):
    """
    Accesses http://charmm-gui.org and uses the Ligand Reader.

    Parameters
    ----------
    file : str
        Path to the input ligand file.
    timeout : int
        Timeout in seconds.
    find_similar_residues : bool
        Whether to tick the "Find similar residues" checkbox before searching.

    Returns
    -------
    filename_output : str
        The absolute path to the output TGZ archive.
    """
    file = _os.path.abspath(file)
    options = _options.Options()
    options.headless = _PC.HEADLESS_CHARMMGUI

    try:
        driver = _seleniumrequests.Chrome(options=options)
    except _exceptions.WebDriverException:
        try:
            driver = _seleniumrequests.Firefox(options=options)
        except _exceptions.WebDriverException:
            raise SystemError("Need either Chrome or Firefox for CHARMM-GUI "
                              "functionality.")

    _logging.info("Accessing http://www.charmm-gui.org ...")
    driver.get("http://www.charmm-gui.org/?doc=input/ligandrm")

    pdb_element = driver.find_element_by_name("file2")
    pdb_element.send_keys(file)

    upload_button = driver.find_element_by_xpath(
        "//input[@type='button' and @value='Upload MOL/MOL2/SDF']")
    upload_button.click()

    driver.switch_to.alert.accept()

    _time.sleep(5)

    if find_similar_residues:
        checkbox = driver.find_element_by_name("simi")
        checkbox.click()

    next_button = driver.find_element_by_id("nextBtn")
    next_button.click()

    # could add some support for options. For now, we just go with the
    # defaults.
    next_button = driver.find_element_by_id("nextBtn")
    next_button.click()

    try:
        _logging.info("Retrieving files...")
        wait = _wait.WebDriverWait(driver, timeout)
        wait.until(_EC.visibility_of_any_elements_located((_by.By.CLASS_NAME,
                                                           "download")))
    except TimeoutError:
        raise ConnectionError("Could not retrieve any files. Please increase "
                              "the maximum timeout or try again later.")

    _logging.info("Downloading TGZ archive...")
    filebase = _os.path.splitext(file)[0]
    tgz_file = driver.find_elements_by_partial_link_text(".tgz")[0]
    response = driver.request('POST', tgz_file.get_attribute("href"),
                              verify=False, stream=True)
    with open(filebase + "_CHARMM.tgz", "wb") as file:
        file.write(response.raw.read())

    driver.quit()
    return filebase + "_CHARMM.tgz"


def charmmguiTransform(filename, **kwargs):
    """
    Adds missing residues to a PDB file.

    Parameters
    ----------
    filename : str
        Name of the input PDB file.
    kwargs
        Keyword arguments to be passed to PDBReader.

    Returns
    -------
    filename_output : str
        Absolute path to the modified file.
    """
    tgz = _tarfile.open(PDBReader(filename, **kwargs), "r:gz")
    paths = []
    for member in tgz:
        if "model.pdb" in member.name.split("/")[-1]:
            f = tgz.extractfile(member)
            with open(member.name.split("/")[-1], "wb") as file:
                file.write(f.read())
            paths += [_os.path.abspath(member.name.split("/")[-1])]
    tgz.close()
    paths = sorted(paths)

    pdb_obj = _PDB(paths[0])
    chainID = paths[0].split("pro")[-1][0].upper()
    for chain in pdb_obj:
        chain.chainID = chainID

    for path in paths[1:]:
        pdb_add = _PDB(path)
        chainID = paths[0].split("pro")[-1][0].upper()
        for chain in pdb_add:
            chain.chainID = chainID
        pdb_obj += pdb_add

    filename_output = _os.path.splitext(filename)[0] + "_charmmgui.pdb"
    return fixCharmmguiPDB(pdb_obj, filename, filename_output=filename_output)


def fixCharmmguiPDB(pdb_modified, filename_original, filename_output=None):
    """
    Used to regenerate some data lost by CHARMM-GUI.

    Parameters
    ----------
    pdb_modified : ProtoCaller.IO.PDB.PDB
        PDB object containing all fixed residues created by CHARMM-GUI.
    filename_original : str
        Name of the original PDB file.
    filename_output : str
        Name of the fixed output PDB file.

    Returns
    -------
    filename_output : str
        The absolute path to the fixed output PDB file.
    """

    pdb_original = _PDB(filename_original)

    for miss_res in pdb_original.missing_residues:
        filter = "chainID=='{}'&resSeq=={}&iCode=='{}'".format(
            miss_res.chainID, miss_res.resSeq, miss_res.iCode)
        fixed_res = _copy.copy(pdb_modified.filter(filter)[0])
        fixed_res.chainID = miss_res.chainID
        fixed_res.resSeq = miss_res.resSeq
        if fixed_res.resName != miss_res.resName:
            _warnings.warn("Mismatch between original residue name ({}) "
                           "and the residue name output by CHARMM-GUI ({}) in "
                           "chain {}, residue number {}.".format(
                miss_res.resName, fixed_res.resName,
                miss_res.chainID, miss_res.resSeq))
        chain = pdb_original.filter("chainID=='{}'".format(miss_res.chainID),
                                    type="chains")[0]

        if miss_res > chain[-1]:
            chain.append(fixed_res)
        else:
            for i, res in enumerate(chain):
                if res > miss_res:
                    chain.insert(i, fixed_res)
                    break
    pdb_original.missing_residues = []

    pdb_original.reNumberAtoms()
    pdb_original.reNumberResidues()
    if filename_output is None:
        filename_output = _os.path.splitext(pdb_original.filename)[0] + \
                          "_modified.pdb"

    return pdb_original.writePDB(filename_output)
