"""Run an admitted timing window with user-authorized temporary Warp suspension."""
import argparse,json,os,signal,subprocess,sys,time
from pathlib import Path
from campaign import ROOT,EVIDENCE,DOCS,sha,save

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output',required=True)
p.add_argument('--start-at',default='cli_ab')
p.add_argument('--suspend-warp',action='store_true')
a=p.parse_args()
directory=EVIDENCE/(a.output+'_host_control');directory.mkdir(exist_ok=False)
def identity(pid):
    try:
        f=(Path('/proc')/str(pid)/'stat').read_text();parts=f[f.rindex(')')+2:].split()
        return parts[19],parts[0]
    except (OSError,ValueError):return None
receipt={'authorization':'User explicitly resumed benchmarks and permitted freezing Warp if needed, followed by restoration and completion of the campaign.',
         'controller_sha256':sha(Path(__file__)),'suspended':[], 'restored':[],
         'kernel':os.uname().release,'start_at':a.start_at,'complete':False}
save(directory/'receipt.json',receipt)
watchdog_code=r'''
import json,os,signal,sys,time
from pathlib import Path
parent,start,record=sys.argv[1:]
def ident(pid):
    try:
        s=Path('/proc',str(pid),'stat').read_text();return s[s.rindex(')')+2:].split()[19]
    except (OSError,ValueError):return None
deadline=time.monotonic()+10800
while ident(parent)==start and time.monotonic()<deadline:time.sleep(1)
d=json.loads(Path(record).read_text())
for r in d['suspended']:
    if ident(r['pid'])==r['starttime']:
        try:os.kill(r['pid'],signal.SIGCONT)
        except ProcessLookupError:pass
Path(record).with_name('watchdog_restore.json').write_text(json.dumps({'restored_after_parent_exit_or_deadline':True,'processes':d['suspended']}))
'''
watchdog=None;child=None
def interrupted(signum,frame):
    if child is not None and child.poll() is None:child.send_signal(signal.SIGTERM)
    raise SystemExit(128+signum)
for sig in [signal.SIGTERM,signal.SIGINT,signal.SIGHUP]:signal.signal(sig,interrupted)
try:
    if a.suspend_warp:
        candidates=[]
        for process in Path('/proc').iterdir():
            if not process.name.isdigit():continue
            try:
                if (process/'comm').read_text().strip()!='warp-terminal':continue
                ident=identity(int(process.name))
                if ident and ident[1] not in ['T','t']:candidates.append({'pid':int(process.name),'starttime':ident[0]})
            except OSError:continue
        receipt['suspended']=candidates;save(directory/'receipt.json',receipt)
        # Independent session restores the GUI even if this controller is killed.
        watchdog=subprocess.Popen([sys.executable,'-c',watchdog_code,str(os.getpid()),identity(os.getpid())[0],str(directory/'receipt.json')],start_new_session=True,stdin=subprocess.DEVNULL,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        receipt['watchdog_pid']=watchdog.pid;save(directory/'receipt.json',receipt)
        for row in candidates:
            if identity(row['pid']) and identity(row['pid'])[0]==row['starttime']:os.kill(row['pid'],signal.SIGSTOP)
        receipt['states_after_suspend']={str(r['pid']):identity(r['pid']) for r in candidates};save(directory/'receipt.json',receipt)
        time.sleep(15)
    cmd=[sys.executable,str(DOCS/'performance_window.py'),'--output',a.output,'--start-at',a.start_at]
    receipt['command']=cmd;save(directory/'receipt.json',receipt)
    child=subprocess.Popen(cmd,cwd=ROOT)
    receipt['returncode']=child.wait();receipt['complete']=True
finally:
    for row in receipt['suspended']:
        ident=identity(row['pid'])
        if ident and ident[0]==row['starttime']:
            os.kill(row['pid'],signal.SIGCONT);receipt['restored'].append(row)
    save(directory/'receipt.json',receipt)
    if watchdog is not None and watchdog.poll() is None:
        watchdog.terminate();watchdog.wait(timeout=5)
raise SystemExit(receipt.get('returncode',1))
