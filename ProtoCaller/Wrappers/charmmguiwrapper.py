import os as _os
import tarfile as _tarfile
import time as _time

from selenium.common import exceptions as _exceptions
from selenium.webdriver.firefox import options as _options
from selenium.webdriver.common import by as _by
from selenium.webdriver.support import expected_conditions as _EC
from selenium.webdriver.support import wait as _wait
from ProtoCaller.shared import seleniumrequests as _seleniumrequests

from ProtoCaller.IO.PDB import PDB as _PDB


def PDBReader(file, timeout=60):
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
    file = _os.path.abspath(file)
    options = _options.Options()
    options.headless = True

    try:
        driver = _seleniumrequests.Chrome(options=options)
    except _exceptions.WebDriverException:
        try:
            driver = _seleniumrequests.Firefox(options=options)
        except _exceptions.WebDriverException:
            raise SystemError("Need either Chrome or Firefox for CHARMM-GUI "
                              "functionality.")

    print("Accessing http://www.charmm-gui.org ...")
    driver.get("http://www.charmm-gui.org/?doc=input/pdbreader")

    pdb_element = driver.find_element_by_name("file")
    pdb_element.send_keys(file)

    pdb_radio = driver.find_element_by_xpath("//input[@name='pdb_format' and "
                                             "@value='PDB']")
    pdb_radio.click()

    next_button = driver.find_element_by_id("nextBtn")
    next_button.click()

    # could add some support for options. For now, we just go with the
    # defaults.
    next_button = driver.find_element_by_id("nextBtn")
    next_button.click()

    next_button = driver.find_element_by_id("nextBtn")
    next_button.click()

    try:
        print("Retrieving files...")
        wait = _wait.WebDriverWait(driver, timeout)
        wait.until(_EC.visibility_of_any_elements_located((_by.By.CLASS_NAME,
                                                           "download")))
    except TimeoutError:
        raise ConnectionError("Could not retrieve any files. Please increase "
                              "the maximum timeout or try again later.")

    print("Downloading TGZ archive...")
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
    options.headless = True

    try:
        driver = _seleniumrequests.Chrome(options=options)
    except _exceptions.WebDriverException:
        try:
            driver = _seleniumrequests.Firefox(options=options)
        except _exceptions.WebDriverException:
            raise SystemError("Need either Chrome or Firefox for CHARMM-GUI "
                              "functionality.")

    print("Accessing http://www.charmm-gui.org ...")
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
        print("Retrieving files...")
        wait = _wait.WebDriverWait(driver, timeout)
        wait.until(_EC.visibility_of_any_elements_located((_by.By.CLASS_NAME,
                                                           "download")))
    except TimeoutError:
        raise ConnectionError("Could not retrieve any files. Please increase "
                              "the maximum timeout or try again later.")

    print("Downloading TGZ archive...")
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

    for path in paths[1:]:
        pdb_obj += _PDB(path)

    return fixCharmmguiPDB(pdb_obj, filename)


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

    for missing_residue in pdb_original.missing_residues:
        modelled_res = pdb_modified.filter("chainID=='{}'&resSeq=={}".format(
            missing_residue.chainID, missing_residue.resSeq))[0]

        for i, chain in enumerate(pdb_original):
            breakloops = False
            if chain.type == "chain":
                for j, residue in enumerate(chain):
                    if all([j == 0, modelled_res.chainID == residue.chainID,
                            modelled_res < residue]):
                        chain.insert(j, modelled_res)
                        breakloops = True
                    elif all([j == len(chain) - 1,
                              modelled_res.chainID == residue.chainID,
                              residue < modelled_res]):
                        chain.insert(j + 1, modelled_res)
                        breakloops = True
                    elif all([0 < j < len(chain),
                              chain[j - 1] < modelled_res < chain[j]]):
                        chain.insert(j, modelled_res)
                        breakloops = True
                    if breakloops:
                        break
            if breakloops:
                break
    pdb_original.missing_residues = []

    pdb_original.reNumberAtoms()
    if filename_output is None:
        filename_output = _os.path.splitext(pdb_original.filename)[0] + \
                          "_modified.pdb"
    pdb_original.writePDB(filename_output)

    return _os.path.abspath(filename_output)
