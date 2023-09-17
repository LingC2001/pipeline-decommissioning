
import sys
import socket
import selectors
import types
from utils import *
import random as rnd

FILE_DEST = '../../Databases/camera/'
MESSAGE_DEST = '../../Databases/status_data/'

class ServerConnection:

    def __init__(self, host, port): 
        # TODO: unique identifier
        self.id = rnd.randint(0,1000)
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.host = host
        self.port = port
        self.header = None
        self.mode = None
        self.filename = None
        
class Server:

    def __init__(self, host, ports):
        self.connections = {} # maps sockets to ServerConnections
        self.selector = selectors.DefaultSelector()
        self.host = host
        self.ports = ports

    def open_socket(self, port):
        # Initiate a TCP socket
        new_sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        new_sock.bind((self.host, port))
        new_sock.listen()
        new_sock.setblocking(False)
        return new_sock
    
    def register_socket(self, ports):
        # Open a TCP socket
        for port in ports:
            sock = self.open_socket(port)
            self.selector.register(sock, selectors.EVENT_READ, data=None)
            print(f"Listening on {(self.host, port)}")

    def accept_connection(self, sock):
        conn, addr = sock.accept()
        new_conn = ServerConnection (
            host = addr[0],
            port = addr[1]
        )
        self.connections [conn] = new_conn
          # Should be ready to read
        print(f"Accepted connection from {addr}")
        conn.setblocking(False)
        data = types.SimpleNamespace(
            addr=addr,
            inb=b"",
            outb=b""
            )
        events = selectors.EVENT_READ | selectors.EVENT_WRITE
        self.selector.register(conn, events, data=data)


    def service_connection(self, key, mask):

        sock = key.fileobj
        data = key.data
        
        if mask & selectors.EVENT_READ:
            recv_data = sock.recv(1024)  # Should be ready to read            

            curr_conn = self.connections[sock]

            if not curr_conn.header:
                print("Received header")
                header, data = recv_data.split(b"$SPLIT$")
                curr_conn.header = eval(header)
                curr_conn.mode = curr_conn.header['mode']
                curr_conn.filename = curr_conn.header['filename']
                recv_data = data

            if recv_data and curr_conn.header:
                
                if curr_conn.header['mode']== 'files' :
                    f_name = curr_conn.header['filename'].split('/')[-1]
                    write_file (FILE_DEST + f_name, recv_data)

                elif curr_conn.header['mode'] in ['messages', 'dust']:
                    store_message (MESSAGE_DEST+f"{curr_conn.header['mode']}.csv", recv_data.decode())

                elif 'robot' in curr_conn.header['mode'] :
                    store_message (MESSAGE_DEST+f"{curr_conn.header['mode']}.csv", recv_data.decode())
                
                elif curr_conn.header['mode'] == 'coords' :
                    parsed_data = recv_data.decode()[1:-1]
                    store_message (MESSAGE_DEST+f"{curr_conn.header['mode']}.csv", parsed_data.decode())

                elif curr_conn.header['mode'] == 'IR' :
                    parsed_data = recv_data.decode()[1:-1]
                    store_message (MESSAGE_DEST+f"{curr_conn.header['mode']}.csv", parsed_data)
                else:
                    raise TypeError ('The 4th keword arg. "mode", must be file or message.')

            if not recv_data : # empty message found in the buffer
                self.selector.unregister(sock)
                sock.close()

    def run(self):
        # TODO: Enable listening to multiple ports
        # TODO : timeout connection
        # print("Starting the server ...")
        # print("-------------------------------------------------------")
        # inp_ = input ("Please specify a list of listening ports: (0-65536)\n")
        # ports = eval(inp_)
        # print("----------")
        self.register_socket(self.ports)
        try:
            while True:
                events = self.selector.select(timeout=None)
                for key, mask in events:
                    if key.data is None: 
                        self.accept_connection(key.fileobj) # Start a new connection
                    else:
                        self.service_connection(key, mask)
        except KeyboardInterrupt:
            print("Caught keyboard interrupt, exiting")
        finally:
            self.selector.close()
            



if __name__ == '__main__' :
    #Start the client
    ports = eval(sys.argv[1])
    server = Server('127.0.0.1', ports)
    server.run()
