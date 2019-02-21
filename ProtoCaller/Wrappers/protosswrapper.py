import os as _os
import re as _re
import warnings as _warnings

from selenium.webdriver.firefox import options as _options
from selenium.webdriver.common import by as _by
from selenium.webdriver.support import expected_conditions as _EC
from selenium.webdriver.support import wait as _wait
from ProtoCaller.shared import seleniumrequests as _seleniumrequests


def protossTransform(filename_pdb, filename_sdf=None, timeout=60, relative=True):
    if relative == True:
        filename_pdb = _os.getcwd() + "/" + filename_pdb
        if filename_sdf is not None: filename_sdf = _os.getcwd() + "/" + filename_sdf

    options = _options.Options()
    options.set_headless(headless=True)

    try:
        driver = _seleniumrequests.Chrome(chrome_options=options)
    except:
        try:
            driver = _seleniumrequests.Firefox(firefox_options=options)
        except:
            print("Need either Chrome or Firefox for Protoss functionality. Otherwise, pathway might be corrupt")
            return -1

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
        wait.until(_EC.visibility_of_any_elements_located((_by.By.ID, "protossdownloadpdb")))
    except TimeoutError:
        print("Could not retrieve any files. Please increase maximum timeout or try later.")
        return -1

    ligand_address, pdb_address, log_address = driver.find_elements_by_css_selector("[action*=download]")
    pdbCode = _re.search(r"proteins.plus/([^/]*)/", pdb_address.get_attribute("action")).group(1)
    post_params = {"pdbCode" : pdbCode}

    print("Downloading files...")
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore")
        filenames = ["protein_protoss.pdb", "ligands_protoss.sdf", "log_protoss.pdb"]
        for link, filename in zip([pdb_address, ligand_address, log_address], filenames):
            if link is ligand_address and filename_sdf is None:
                filenames[1] = None
                continue
            response = driver.request('POST', link.get_attribute("action"), data=post_params, verify=False)
            with open(filename, "w") as file:
                for line in response.content.decode():
                    file.write(line)

    driver.quit()
    #fix annoying numbering issue with protoss
    filenames[0] = fixProtossPDB(filenames[0], filenames[0])
    return filenames

def fixProtossPDB(filename_input, filename_output=None):
    if filename_output is None:
        filename_output = filename_input[:-4] + "_modified.pdb"

    file_input = open(filename_input).readlines()
    with open(filename_output, "w") as file_output:
        for line in file_input:
            if line[:6] == "HETATM" and (line[25].isalpha() or line[25] == " "):
                line = line[:22] + " " + line[22:26] + line[27:]
            file_output.write(line)

    return filename_output
