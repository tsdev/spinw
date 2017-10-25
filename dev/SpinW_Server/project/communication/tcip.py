import socket

class tcip:
    def __init__(self, host, port=13001):
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.port = port
        self.host = host
        self.is_connected = False

    def __del__(self):
        if self.is_connected:
            self.socket.close()

    def connect(self):
        if not self.is_connected:
            self.socket.connect((self.host, self.port))
            self.is_connected = True

    def send_comand(self, cmd):
        try:
            if self.is_connected:
                self.socket.sendall(cmd.encode('utf-8'))
            else:
                self.connect()
                self.socket.sendall(cmd.encode('utf-8'))
            return True
        except:
            self.socket.connect((self.host, self.port))
            self.connect()
            self.socket.sendall(cmd.encode('utf-8'))
            return False