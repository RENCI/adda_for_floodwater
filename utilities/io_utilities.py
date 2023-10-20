#!/usr/bin/env python

import sys,os
import yaml
import json

def construct_base_rootdir(inconfig=None, base_dir_extra='None')->str:
    """
    Report the TOP of the desired data tree usually specified as $RUNTIMEDIR
    The data is generally read from the main.yml
    Create the directory tree is needed

    Parameters:
        inconfig str: usually read from main.yml
            Must either result from  config['DEFAULT']['RDIR'] or be a valid str to a local dir ,
        base_dir_extra <str>. A convenience feature that allow a user to specify a fixed output tree
            specified in $RUNTIMEDIR or a specific dir but then build (for example) individuakl experimental nametags underneath
    Returns:
        rootdir <str>: as rootdir
    """
    try:
        rootdir = os.environ[inconfig.replace('$', '')]  # Yaml call to be subsequently removed
    except OSError:
       print('Chosen basedir invalid: check '+str(inconfig['DEFAULT']['RDIR']))
    except KeyError:
       rootdir=inconfig
    if base_dir_extra is not None:
        rootdir = rootdir+'/'+base_dir_extra
    if not os.path.exists(rootdir):
        try:
            os.makedirs(rootdir)
        except OSError:
            print("Creation of the high level run directory %s failed" % rootdir)
    return rootdir

def get_full_filename_with_subdirectory_prepended(rootdir, new_subdir, fname )->str:
    """
    It is expected that rootdir already exists
    subdir is created as desired

    Parameters:
        rootdir <str>: usually the current TOP level  as returned by contruct_base_rootdir
        new_subdir <str>: Optional new appended subdir
        fname <str>: filename with extension
    Returns: 
        rootdir/new_subdir/fname full path for subsequent write by the caller
    """

    if not os.path.exists(rootdir):
        print('rootdir must exist. It doesnt: Try calling contruct_base_rootdir(inconfig) first. %s' % rootdir)
        sys.exit(1)
    fulldir = os.path.join(rootdir, new_subdir)
    if not os.path.exists(fulldir):
        try:
            os.makedirs(fulldir)
        except OSError:
            if not os.path.isdir(fulldir): 
                print("Creation of the rootdir/new_subdir directory %s failed" % fulldir)
                sys.exit(1)
    return os.path.join(fulldir, fname)

def write_pickle(df, rootdir=None,subdir=None,fileroot=None,iometadata='')->str:
    """ 
    Write the indicated dataframe to the file
    rootdir/subdir/fileroot+iometadata.pkl

    Parameters:
        df <DataFrame>: of data to write to disk
        rootdir <str>: top of the output FS tree
        subdir <str>: additional subdirectory into which to place file
        fileroot <str : filename with no extension
        iometadata <str>: additional data to append to fileroot before extension.

    Returns full filename for written data
    """
    fileroot = '_'.join([fileroot,iometadata]) if iometadata != '' else fileroot
    try:
        newfilename = get_full_filename_with_subdirectory_prepended(rootdir, subdir, fileroot+'.pkl')
        df.to_pickle(newfilename)
        #print('Wrote PICKLE file {}'.format(newfilename))
    except OSError:
        raise OSError("Failed to write PKL file %s" % (newfilename))
    return newfilename

def write_csv(df, rootdir=None,subdir=None,fileroot=None,iometadata='')->str:
    """ 
    Write the indicated dataframe to the file
    rootdir/subdir/fileroot+iometadata.csv

    Parameters:
        df <DataFrame>: of data to write to disk
        rootdir <str>: top of the output FS tree
        subdir <str>: additional subdirectory into which to place file
        fileroot <str : filename with no extension
        iometadata <str>: additional data to append to fileroot before extension.

    Returns full filename for written data
    """
    fileroot = '_'.join([fileroot,iometadata]) if iometadata != '' else fileroot
    try:
        newfilename = get_full_filename_with_subdirectory_prepended(rootdir, subdir, fileroot+'.csv')
        df.to_csv(newfilename)
        #print('Wrote CSV file {}'.format(newfilename))
    except OSError:
        raise OSError("Failed to write CSV file %s" % (newfilename))
    return newfilename


def write_json(df, rootdir=None,subdir=None,fileroot=None,iometadata='')->str:
    """ 
    Write the indicated dataframe to the file
    rootdir/subdir/fileroot+iometadata.json

    Parameters:
        df <DataFrame>: of data to write to disk
        rootdir <str>: top of the output FS tree
        subdir <str>: additional subdirectory into which to place file
        fileroot <str : filename with no extension
        iometadata <str>: additional data to append to fileroot before extension.

    Returns full filename for written data
    """
    fileroot = '_'.join([fileroot,iometadata]) if iometadata != '' else fileroot
    try:
        newfilename = get_full_filename_with_subdirectory_prepended(rootdir, subdir, fileroot+'.json')
        df.to_json(newfilename)
        #print('Wrote JSON file {}'.format(newfilename))
    except OSError:
        raise OSError("Failed to write JSON file %s" % (newfilename))
    return newfilename


def write_dict_to_json(dictdata, rootdir=None,subdir=None,fileroot=None,iometadata='')->str:
    """ 
    Write the indicated dataframe to the file
    rootdir/subdir/fileroot+iometadata.json

    If the index/keys are datetime (or equiv) you must have converted thenm to strings

    Parameters:
        dictdata <dict>: of data to write to disk
        rootdir <str>: top of the output FS tree
        subdir <str>: additional subdirectory into which to place file
        fileroot <str : filename with no extension
        iometadata <str>: additional data to append to fileroot before extension.

    Returns full filename for written data
    """
    fileroot = '_'.join([fileroot,iometadata]) if iometadata != '' else fileroot
    try:
        newfilename = get_full_filename_with_subdirectory_prepended(rootdir, subdir, fileroot+'.json')
        with open(newfilename, 'w') as fp:
            json.dump(dictdata, fp)
        #print('Wrote dictdata to JSON file {}'.format(newfilename))
    except OSError:
        raise OSError("Failed to write dictdata to JSON file %s" % (newfilename))
    return newfilename

def read_json_file(filepath)->dict:
    # Read data from JSON file specified by full path
    data = {}
    try:
        with open(filepath, 'r') as fp:
            data = json.load(fp)
    except FileNotFoundError:
        raise FileNotFoundError("Failed to read file %s" % (filepath))
    return data

def write_json_file(data, filepath):
    # write data from JSON file specified by full path
    try:
        with open(filepath, 'w') as fp:
            json.dump(data, fp)
    except OSError:
        raise OSError("Failed to write JSON file %s" % (filepath))

# Special case output that is often used in our EDS context

def write_ADCIRC_formatted_gridfield_to_Disk(df, value_name='VAL',rootdir=None,subdir=None,filename=None,iometadata='')->str:
    """
    Write that generally is an adcirc irregular grid with values. Generally, this might be an extrapolated offset field.
    The output file is suitable for reading by ADCIRC.

    Parameters:
        dfx <DataFrame> with columns: (lon,lat,value).

    Results:
        df  saved to rootdir/interpolated/ADCIRC_interpolated_wl_metadata.csv
    """
    #fileroot = '_'.join([fileroot,iometadata]) if iometadata != '' else fileroot
    #newfilename = get_full_filename_with_subdirectory_prepended(rootdir, subdir, fileroot)
    newfilename = filename
    df_adcirc = df[value_name].to_frame().astype(str)
    df_adcirc['node']=(df_adcirc.index+1).astype(str) # NODEID is index id +1
    d = []
    d.append('# Interpolated field')
    d.append('99999.9')
    d.append('0.0')
    for index,row in df_adcirc.iterrows():
        nd = row['node']
        nv = row[value_name]
        d.append(nd+','+nv)
    with open(newfilename, mode='wt', encoding='utf-8') as myfile:
        myfile.write('\n'.join(d))
    #print('Wrote current extrapolated ADCIRC grid to disk')
    return newfilename
