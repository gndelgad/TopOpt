import fileinput
import sys
import re

def replace(file, pattern, subst):
    # Read contents from file as a single string
    file_handle = open(file, 'rb')
    file_string = file_handle.read()
    file_handle.close()

    # Use RE package to allow for replacement (also allowing for (multiline) REGEX)
    file_string = (re.sub(pattern, subst, file_string))

    # Write contents to file.
    # Using mode 'w' truncates the file.
    file_handle = open(file, 'wb')
    file_handle.write(file_string)
    file_handle.close()
	
replace("Th.LevelSet.mesh","MeshVersionFormatted 0","MeshVersionFormatted 1")