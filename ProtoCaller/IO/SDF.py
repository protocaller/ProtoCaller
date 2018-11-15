def mergeSDFs(input_files, output_filename=None):
    if input_files in [None, []]:
        return None
    if output_filename is None: output_filename = "ligands_merged.sdf"

    with open(output_filename, "w") as output:
        for file in input_files:
            content = open(file).readlines()
            for line in content:
                output.write(line)

    return output_filename

def splitSDFs(input_files):
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
