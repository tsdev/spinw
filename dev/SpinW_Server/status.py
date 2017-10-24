import threading
from project import db, config, results
from project.models import User, UserJobs
import datetime
import os
import time

class Status(threading.Thread):
    def __init__(self, delay):
        threading.Thread.__init__(self)
        self.delay = delay

    def run(self):
        while True:
            jobs = db.session.query(UserJobs).filter(UserJobs.completed == False, UserJobs.running == True).order_by(UserJobs.start_time).all()
            for job in jobs:
                token = job.token
                if config.get('USE_PYMATLAB'):
                    if token in results:
                        isR = results[token].done()
                    else:
                        isR = False
                else:
                    if os.path.isfile(os.path.join(config.get('UPLOAD_FOLDER'), 'out_%s.mat' % token)):
                        isR = True
                    else:
                        isR = False
                if isR:
                    job.completed = True
                    job.end_time = datetime.datetime.now()
                    user = db.session.query(User).filter_by(id=job.user_id).first()
                    user.quota_used = user.quota_used + (job.start_time - job.end_time).total_seconds()
                    db.session.commit()
            time.sleep(self.delay)