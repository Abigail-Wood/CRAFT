import sys

from .main import main

# A hack to enable python -m craft --more --switches --go --here
# lifted from unittest/_main_.py with tweaks.
if sys.argv[0].endswith("__main__.py"):
    import os.path
    # change sys.argv[0] to make help message more useful.
    # Use executable without path, unquoted
    # (it's just a hint anyway)
    # If you have spaces in your executable you get what you deserve!
    executable = os.path.basename(sys.executable)
    sys.argv[0] = executable + " -m " + __spec__.parent
    del os

if __name__=='__main__':
	main()
