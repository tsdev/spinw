#################
#### imports ####
#################

import datetime
from project import config, db
from passlib.apps import custom_app_context as pwd_context
from itsdangerous import (TimedJSONWebSignatureSerializer
                          as Serializer, BadSignature, SignatureExpired)

if config.get('USE_LDAP'):
    import ldap
    ldap.set_option(ldap.OPT_REFERRALS, 0)
    ldap.protocol_version = 3
    conn = ldap.initialize('ldap://%s' % config.get('LDAP_SERVER'))
else:
    conn = []

class AutoSerialize(object):
    'Mixin for retrieving public fields of model in json-compatible format'
    __public__ = None

    def get_public(self, exclude=(), extra=()):
        "Returns model's PUBLIC data for jsonify"
        data = {}
        keys = self._sa_instance_state.attrs.items()
        public = self.__public__ + extra if self.__public__ else extra
        for k, field in keys:
            if public and k not in public: continue
            if k in exclude: continue
            value = self._serialize(field.value)
            if value is not None:
                data[k] = value
            else:
                data[k] = 0
        return data

    @classmethod
    def _serialize(cls, value, follow_fk=False):
        if type(value) is datetime.datetime:
            ret = value.strftime("%Y-%m-%d %H:%M:%S")
        elif isinstance(value, str):
            ret = value
        elif hasattr(value, '__iter__'):
            ret = []
            for v in value:
                ret.append(cls._serialize(v))
        elif AutoSerialize in value.__class__.__bases__:
            ret = value.get_public()
        else:
            ret = value

        return ret


class User(db.Model, AutoSerialize):
    __tablename__ = 'users'
    __public__ = ('id', 'username', 'registered_on', 'admin', 'quota_total', 'quota_used')
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(32), index=True, unique=True)
    token = db.Column(db.String(10))
    jobs = db.relation('UserJobs', backref='users', lazy='dynamic')
    password_hash = db.Column(db.String(64))
    registered_on = db.Column(db.DateTime, nullable=False)
    admin = db.Column(db.Boolean, nullable=False, default=False)
    confirmed = db.Column(db.Boolean, nullable=False, default=False)
    confirmed_on = db.Column(db.DateTime, nullable=True)
    quota_total = db.Column(db.Float)
    quota_used = db.Column(db.Float)

    def __init__(self, username, password, useLDAP=False, admin=False, confirmed_on=None, confirmed=False):
        self.username = username
        if useLDAP:
            self.password_hash = None
        else:
            self.hash_password(password)
        self.registered_on = datetime.datetime.now()
        self.admin = admin
        self.confirmed = confirmed
        self.confirmed_on = confirmed_on
        self.quota_total = float("inf")
        self.quota_used = 0
        self.jobs = []
    def hash_password(self, password):
        self.password_hash = pwd_context.encrypt(password)

    def verify_password(self, password):
        try:
            if config.get('useLDAP'):
                conn.simple_bind_s('%s@%s' % (self.username, config.get('LDAP_SERVER')), password)
                return True
            else:
                return pwd_context.verify(password, self.password_hash)
        except:
            return False

    def generate_auth_token(self, expiration=600):
        s = Serializer(config['SECRET_KEY'], expires_in=expiration)
        return s.dumps({'id': self.id}, salt=config['SECURITY_PASSWORD_SALT'])
        # return s.dumps({'id': self.id})

    @staticmethod
    def verify_auth_token(token):
        s = Serializer(config['SECRET_KEY'])
        try:
            data = s.loads(token, salt=config['SECURITY_PASSWORD_SALT'])
            # data = s.loads(token)
        except SignatureExpired:
            return None  # valid token, but expired
        except BadSignature:
            return None  # invalid token
        user = User.query.get(data['id'])
        return user


class UserJobs(db.Model, AutoSerialize):
    __tablename__ = 'jobs'
    __public__ = ('job_id', 'token', 'start_time', 'end_time', 'completed')
    job_id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey("users.id"))
    token = db.Column(db.String(32), index=True)
    start_time = db.Column(db.DateTime)
    end_time = db.Column(db.DateTime)
    running = db.Column(db.Boolean)
    completed = db.Column(db.Boolean)

    def __init__(self, token):
        self.token = token
        self.completed = False
        self.running = False