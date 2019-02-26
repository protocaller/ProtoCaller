def mergeSDFs(input_files, output_filename="ligands_merged.sdf"):
    """
    Merges all input SDF files into a single one.

    Parameters
    ----------
    input_files : [str]
        Names of the input SDF files.
    output_filename: str
        Name of the output file.

    Returns
    -------
    output_filename: str
        Absolute path to the output file.
    """
    if input_files in [None, []]:
        return None

    with open(output_filename, "w") as output:
        for file in input_files:
            content = open(file).readlines()
            for line in content:
                output.write(line)

    return output_filename

def splitSDFs(input_files):
    """
    Splits each of the input SDF files into many SDF files, each containing only one molecule.

    Parameters
    ----------
    input_files : [str]
        Names of the input SDF files.

    Returns
    -------
    output_files: [str]
        Absolute paths to the output files.
    """
    if input_files in [[], None, [None]]:
        return []
    output_files = []
    for input_file in input_files:
        content = open(input_file).readlines()
        beginning = True
        for line in content:
            if beginning:
                output_files += [line.strip() + "_split.sdf"]
                beginning = False
                output_file = open(output_files[-1], "w")
            if line.strip() == "$$$$":
                beginning = True
            output_file.write(line)
    return output_files
