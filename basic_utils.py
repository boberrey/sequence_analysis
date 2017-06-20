#!/usr/bin/env python

"""
Basic utility functions for sequence analysis scripts

Ben Ober-Reynolds
"""

import os

def find_files_in_directory(dirPath, extensionList=None, 
                            excludedExtensionList=None):
    """
    Locate files in a given directory path. Optionally, desired files are 
    identified as matching one of the extension types provided in 
    'extensionList'
    Input: directory path, list of approved extensions, (list of excluded 
        extensions)
    Output: List of found files 
    """
    def extension_match(filename, extensionList=None):
        # from CPlibs
        if extensionList is not None:
            for currExt in extensionList:
                if filename.lower().endswith(currExt.lower()):
                    return True
        return False

    dirList = os.listdir(dirPath)
    fileList = []
    for currFilename in dirList:
        if (extension_match(currFilename, extensionList) 
        and not extension_match(currFilename, excludedExtensionList)): 
            fileList.append(dirPath+currFilename)
    if len(dirList) == 0:
        print('\tNONE FOUND')
    else:
        return fileList