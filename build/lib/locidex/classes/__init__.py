from subprocess import Popen, PIPE
def run_command(command):
    """
    This method runs the passed command on the shell command line.

    Parameters
    ----------
    command : str
        The command for the command line

    Returns
    -------
    bytes, bytes
        stdout, stderr: the standard out and error messages returned by the command line
    """
    p = Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    stdout = stdout.decode('utf-8')
    stderr = stderr.decode('utf-8')
    return stdout, stderr
