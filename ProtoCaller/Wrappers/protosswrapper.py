import os as _os
import re as _re
import warnings as _warnings

from selenium.common import exceptions as _exceptions
from selenium.webdriver.firefox import options as _options
from selenium.webdriver.common import by as _by
from selenium.webdriver.support import expected_conditions as _EC
from selenium.webdriver.support import wait as _wait
from ProtoCaller.shared import seleniumrequests as _seleniumrequests

from ProtoCaller.IO.PDB import PDB as _PDB

__all__ = ["protossTransform"]


def protossTransform(filename_pdb, filename_sdf=None, timeout=60):
    """
    Protonates a protein complex from a PDB and an SDF files using ProToss.

    Parameters
    ----------
    filename_pdb : str
        Name of input PDB file.
    filename_sdf : str
        Name of input SDF file.
    timeout : float
        Timeout in seconds.

    Returns
    -------
    filenames: [str]
        Names of the protonated files.
    """
    filename_pdb = _os.path.abspath(filename_pdb)
    if filename_sdf is not None: filename_sdf = _os.path.abspath(filename_sdf)

    options = _options.Options()
    options.headless = True

    try:
        driver = _seleniumrequests.Chrome(options=options)
    except _exceptions.WebDriverException:
        try:
            driver = _seleniumrequests.Firefox(options=options)
        except _exceptions.WebDriverException:
            raise SystemError("Need either Chrome or Firefox for Protoss "
                              "functionality.")

    print("Accessing https://proteins.plus/ ...")
    driver.get("https://proteins.plus/")

    pdb_element = driver.find_element_by_id("pdb_file_pathvar")
    pdb_element.send_keys(filename_pdb)

    if filename_sdf is not None:
        sdf_element = driver.find_element_by_id("pdb_file_userligand")
        sdf_element.send_keys(filename_sdf)

    go_button = driver.find_element_by_name("commit")
    go_button.click()

    protoss_button = driver.find_elements_by_css_selector("[href*=protoss]")[2]
    protoss_button.click()

    calculate_button = driver.find_element_by_name("commit")
    calculate_button.click()

    try:
        print("Retrieving download links...")
        wait = _wait.WebDriverWait(driver, timeout)
        wait.until(_EC.visibility_of_any_elements_located(
            (_by.By.ID, "protossdownloadpdb")))
    except TimeoutError:
        raise ConnectionError("Could not retrieve any files. Please increase "
                              "the maximum timeout or try again later.")

    sdf, pdb = driver.find_elements_by_css_selector("[action*=download]")[:2]
    pdbCode = _re.search(r"proteins.plus/([^/]*)/",
                         pdb.get_attribute("action")).group(1)
    post_params = {"pdbCode" : pdbCode}

    print("Downloading files...")
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore")

        filenames = [_os.path.splitext(filename_pdb)[0] + "_protoss.pdb"]
        if filename_sdf is None:
            filenames += [None]
        else:
            filenames += [_os.path.splitext(filename_sdf)[0] + "_protoss.sdf"]

        for link, filename in zip([pdb, sdf], filenames):
            if filename is not None:
                response = driver.request('POST', link.get_attribute("action"),
                                          data=post_params, verify=False)
                with open(filename, "w") as file:
                    for line in response.content.decode():
                        file.write(line)

    driver.quit()
    filenames[0] = fixProtossPDB(filenames[0], filename_pdb, filenames[0])

    return filenames


def fixProtossPDB(filename_modified, filename_original, filename_output=None):
    """
    Used to regenerate some data lost by Protoss.

    Parameters
    ----------
    filename_modified : str
        Name of the modified PDB file.
    filename_original : str
        Name of the original PDB file.
    filename_output : str
        Name of the fixed output PDB file.

    Returns
    -------
    filename_output : str
        The absolute path to the fixed output PDB file.
    """
    # fix a misaligned numbering issue with Protoss
    file_input = open(filename_modified).readlines()
    with open(filename_modified, "w") as file_output:
        for line in file_input:
            if (line[:6] == "HETATM" or line[:4] == "ATOM") and \
                    (line[25].isalpha() or line[25] == " "):
                line = line[:22] + " " + line[22:26] + line[27:]
            file_output.write(line)

    pdb_original = _PDB(filename_original)
    pdb_modified = _PDB(filename_modified)

    for res_orig in pdb_original.filter("type=='amino_acid'"):
        filter = "chainID=='{}'&resSeq=={}&iCode=='{}'".format(
            res_orig.chainID, res_orig.resSeq, res_orig.iCode)
        res_mod = pdb_modified.filter(filter)[0]
        res_orig.clear()
        res_orig.__init__(res_mod)

    pdb_original.reNumberAtoms()
    if filename_output is None:
        filename_output = _os.path.splitext(pdb_original.filename)[0] + \
                          "_modified.pdb"

    return pdb_original.writePDB(filename_output)
