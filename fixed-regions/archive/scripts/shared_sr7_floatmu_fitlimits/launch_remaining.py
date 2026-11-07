from run_observed import *
if __name__=='__main__':
    preferred=['9.0','11.5','6.0','12.0']
    masses=preferred+[f'{n/10:.1f}' for n in range(60,121) if f'{n/10:.1f}' not in preferred+['7.5']]
    atomic(R/'remaining_runner.json',dict(status='running',pid=os.getpid(),workers=7,masses=masses,started_at=datetime.datetime.now().astimezone().isoformat()))
    with concurrent.futures.ThreadPoolExecutor(max_workers=7) as pool:
        finished=list(pool.map(job,masses))
    atomic(R/'remaining_runner.json',dict(status='completed',results=finished,workers=7,finished_at=datetime.datetime.now().astimezone().isoformat()))
