#################
#### imports ####
#################
import json

from flask import abort, request, g, render_template, Blueprint
import datetime

from project import db, config, auth
from project.models import User

################
#### config ####
################

user_blueprint = Blueprint('user', __name__, )


@auth.verify_password
def verify_password(username_or_token, password):
    # first try to authenticate by token
    user = User.verify_auth_token(username_or_token)
    if not user:
        # try to authenticate with username/password
        user = User.query.filter_by(username=username_or_token).first()
        if not user or not user.verify_password(password):
            return False
    g.user = user
    return True


@user_blueprint.route('/users/add', methods=['POST'])
def new_user():
    if 'MATLAB' in request.headers.environ.get('HTTP_USER_AGENT'):
        username = request.form.get('username')
        password = request.form.get('password')
    else:
        username = request.json.get('username')
        password = request.json.get('password')
    if username is None:
        abort(400)
    else:
        if password is None and not config.get('USE_LDAP'):
            abort(400)  # missing arguments
    if User.query.filter_by(username=username).first() is not None:
        abort(409)  # existing user
    if config.get('USE_LDAP'):  # Create the user
        user = User(username=username, password="")
    else:
        user = User(username=username, password=password)
    db.session.add(user)
    db.session.commit()
    return (json.dumps({'username': user.username}), 201)


@user_blueprint.route('/users/confirm/<token>')
def confirm_email(token):
    user = User.verify_auth_token(token)
    if user is None:
        ret = {"message": 'The confirmation link is invalid or has expired.', "success": False}
        return json.dumps(ret)
    if user.confirmed:
        ret = {"message": 'Account already confirmed. Please login.', "success": True}
    else:
        user.confirmed = True
        user.confirmed_on = datetime.datetime.now()
        db.session.commit()
        ret = {"message": 'You have confirmed your account. Thanks!', "success": True}
    return json.dumps(ret)


@user_blueprint.route('/users/token')
@auth.login_required
def get_auth_token():
    token = g.user.generate_auth_token(config.get("TOKEN_DURATION"))
    return json.dumps({'token': token.decode('ascii'), 'duration': config.get("TOKEN_DURATION")})


@user_blueprint.route('/users/quota')
@auth.login_required
def get_quota():
    user = g.user
    return json.dumps({'used': user.quota_used, 'total': user.quota_total})


@user_blueprint.route("/users/jobs")
@auth.login_required
def job_list():
    user = g.user
    j_list = []
    for job in user.jobs.all():
        j_list.append(job.get_public())
    return render_template('table.html', jobs=j_list)
