
# Flexible interface to some kind of compute cluster
# Should be capable of running both on Fostino in Grenoble and 
# on the STScI science net at a minimum. Ideally extensible to
# other types of job managers as well.



def configure_cluster(**kwargs):
    """ Configure the details of the compute cluster to be used, 
    and perhaps store to some config file in the user's home directory?"

    """
    raise NotImplementedError("Not yet")

def save_cluster_config(filename):
    """ Save cluster configuration information to some file
    based on currently active cluster configuration"""
    raise NotImplementedError("Not yet")

def read_cluster_config(filename):
    """ Read cluster configuration information from some file
    and restore to become the active cluster configuration"""
    raise NotImplementedError("Not yet")

