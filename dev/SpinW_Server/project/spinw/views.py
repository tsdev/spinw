#################
#### imports ####
#################
import json

from flask import abort, request, g, Blueprint, url_for, send_from_directory
import datetime
import random
import string
import os
import sys
from io import StringIO
from pathlib import Path

from project import db, config, auth, eng, results, running
from project.models import UserJobs


################
#### config ####
################

spinw_blueprint = Blueprint('spinw', __name__,)
out_pipe = StringIO()
err_pipe = StringIO()

@spinw_blueprint.route('/spinw')
def get_resource():
    text = "Welcome to the SpinW compute server"
    return text


@spinw_blueprint.route('/spinw/version')
def get_version():
    temp = {'Version': None, 'Deployed': False, 'Authentication': 'local'}
    if not config.get('USE_PYMATLAB'):
        temp['Deployed'] = True
    if config.get('USE_LDAP'):
        temp['Authentication'] = 'ldap'
    temp['Version'] = config.get('SERVER_VERSION')
    return json.dumps(temp)

@spinw_blueprint.route('/spinw/status/<string:token>')
@auth.login_required
def get_status(token):
    user = g.user
    if not user:
        abort(400)
    job = user.jobs.filter_by(token=token).first()
    return_obj = {"status": "running",
                  "standard": None,
                  "error": None,
                  "url": None}
    if not job:
        return_obj['status'] = "No job submitted"
        return json.dumps(return_obj)

    if token not in results:
        return_obj['status'] = "Not running"
        return json.dumps(return_obj)

    if job.completed:
        try:
            if config.get('USE_PYMATLAB'):
                ret = results[token].result()
            else:
                ret = 'out_%s' % token
            return_obj["status"] = 'done'
            return_obj["url"] = url_for('spinw.download_file', token='%s.mat' % ret, _external=True)
        except:
            e = sys.exc_info()[0]
            if config.get('USE_PYMATLAB'):
                results[token].cancel()
            return_obj["status"] = str(e)
    return json.dumps(return_obj)


@spinw_blueprint.route("/spinw/status/download/<string:token>")
@auth.login_required
def download_file(token):
    path = os.path.join(Path(__file__).parents[2], token)
    if os.path.isfile(path):
        return send_from_directory(directory=Path(__file__).parents[2], filename=token)
    else:
        path = os.path.join(config.get('UPLOAD_FOLDER'), '%s' % token)
        if os.path.isfile(path):
            return send_from_directory(directory=os.path.join(Path(__file__).parents[2],config.get('UPLOAD_FOLDER')), filename=token)
        else:
            return json.dumps({'status': 'File not found.'})

@spinw_blueprint.route('/spinw/upload/<string:filename>', methods=['GET', 'POST'])
@auth.login_required
def upload(filename):
    if request.method == 'POST':
        file = request.data
        extension = os.path.splitext(filename)[1]
        if extension != ".mat":
            abort(400)
        token = ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(8))
        user = g.user
        job = UserJobs(token)
        user.jobs.append(job)
        db.session.commit()
        f_name = token + extension
        if not os.path.isdir(config['UPLOAD_FOLDER']):
            os.makedirs(config['UPLOAD_FOLDER'])

        with open(os.path.join(config['UPLOAD_FOLDER'], f_name), 'wb') as w:
            w.write(file)
        running[token] = eng.getArrayFromByteStream(file)
        # except:
        #     print(Exception.args)
        return (
            json.dumps({'username': user.username, 'status': url_for('spinw.get_status', token=token, _external=True)}), 201)


@spinw_blueprint.route('/spinw/compute/<string:filename>', methods=['GET'])
@auth.login_required
def compute_deployed(filename):
    return abort(404)


@spinw_blueprint.route('/spinw/spinwave/<string:filename>', methods=['POST'])
@auth.login_required
def compute_spinwave(filename):
    user = g.user
    if request.method == 'POST':
        data = request.data
        extension = os.path.splitext(filename)[1]
        if extension != ".mat":
            abort(400)
        job = user.jobs.filter(UserJobs.completed == False).order_by(UserJobs.start_time)[0]
        token = job.token
        f_name = 'in_' + str(token) + extension
        with open(os.path.join(config['UPLOAD_FOLDER'], f_name), 'wb') as w:
            w.write(data)
        if config.get('USE_PYMATLAB'):
            try:
                cmd = eng.getArrayFromByteStream(data)
            except:
                print(Exception.args)
            os.remove(os.path.join(config['UPLOAD_FOLDER'], f_name))
        if token is None:
            abort(400)
    else:
        abort(400)
    if config.get('USE_PYMATLAB'):
        Q = cmd['Q']
        cmd.pop('Q')
        cmd['toFile'] = str(random.randint(1000, 10000))
        job.start_time = datetime.datetime.now()
        job.running = True
        results[token] = eng.spinwavefast(running[token], Q, cmd, async=True, stdout=out_pipe, stderr=err_pipe)
    else:
        job.start_time = datetime.datetime.now()
        job.running = True
        success = eng.send_comand('EXEC ' + token + ' 1.23:')
        results[token] = success
    db.session.commit()
    return json.dumps(
        {'Calculating': True, 'Errors': False, 'status': url_for('spinw.get_status', token=token, _external=True)})


@spinw_blueprint.route('/spinw/powspec/<string:filename>', methods=['POST'])
@auth.login_required
def compute_powspec(filename):
    user = g.user
    if request.method == 'POST':
        data = request.data
        extension = os.path.splitext(filename)[1]
        if extension != ".mat":
            abort(400)
        job = user.jobs.filter(UserJobs.completed == False).order_by(UserJobs.start_time)[0]
        token = job.token
        f_name = 'in_' + str(token) + extension
        with open(os.path.join(config['UPLOAD_FOLDER'], f_name), 'wb') as w:
            w.write(data)
        if config.get('USE_PYMATLAB'):
            try:
                cmd = eng.getArrayFromByteStream(data)
            except:
                print(Exception.args)
            os.remove(os.path.join(config['UPLOAD_FOLDER'], f_name))
        if token is None:
            abort(400)
    else:
        abort(400)
    # try:
    if config.get('USE_PYMATLAB'):
        Q = cmd['Q']
        cmd.pop('Q')
        cmd['toFile'] = str(random.randint(1000, 10000))
        job.start_time = datetime.datetime.now()
        job.running = True
        results[token] = eng.powspecfast(running[token], Q, cmd, async=True, stdout=out_pipe, stderr=err_pipe)
    else:
        job.start_time = datetime.datetime.now()
        job.running = True
        success = eng.send_comand('EXEC ' + token + ' 1.23:')
        results[token] = success
    db.session.commit()
    return json.dumps(
        {'Calculating': True, 'Errors': False, 'status': url_for('get_status', token=token, _external=True)})
