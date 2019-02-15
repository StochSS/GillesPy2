import sys, os, tempfile, shutil

src_dir = os.path.dirname(__file__)
os.chdir(src_dir)
dest_dir = tempfile.mkdtemp()

for files in os.listdir('.'):
    if files != __file__:
        shutil.copy(files,dest_dir)

os.chdir(dest_dir)
os.system("python3 ./run_tests.py")


