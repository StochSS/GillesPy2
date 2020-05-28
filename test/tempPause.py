# import subprocess
# import time
# import signal
# import os
args = ['python3', 'pause_model.py']
# process = subprocess.Popen(args, stdout=subprocess.PIPE)
# time.sleep(2)
# process.send_signal(signal.SIGINT)
# process.kill()
# out, err = process.communicate()
# print(out)


import subprocess
import signal
import time
import os

p = subprocess.Popen(args, preexec_fn=os.setsid,stdout=subprocess.PIPE)
print('sleeping')
time.sleep(4)
os.kill(p.pid, signal.SIGINT)
print('interrupt')
out, err = p.communicate(timeout=1)
print(out)
print('process finished')

