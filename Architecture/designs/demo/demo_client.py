import sys
import socket
import selectors
import types
import random as rnd
import pickle
import numpy as np

class ClientConnection:

    def __init__(self, id_, host, port, message): 
        # TODO: unique identifier
        self.id = id_
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.host = host
        self.port = port
        self.send_registry = message
        self.outbound_data = None
        self.inbound_data = []


class Client:

    def __init__(self, port):
        self.connections = []
        self.selector = selectors.DefaultSelector()
        self.port = int(port)

    # TODO: Move this into ClientConnection
    def initiate_connection (self, host, port, message_packet) :
        server_addr = (host, port)
        print(f"Starting connection {message_packet.connId} to {server_addr}")
        new_con = ClientConnection (
            id_ = message_packet.connId,
            host = host,
            port = port,
            message = message_packet,
        )
        new_con.socket.setblocking(False)
        new_con.socket.connect_ex(server_addr)
        events = selectors.EVENT_READ | selectors.EVENT_WRITE
        self.connections.append(new_con)
        self.selector.register(new_con.socket, events, data=message_packet)


    def service_connection(self, key, mask):
        sock = key.fileobj
        data = key.data

        if mask & selectors.EVENT_READ:
            recv_data = sock.recv(1024)  # Should be ready to read
            if recv_data:
                print(f"Client : Received {recv_data!r} from connection {data.connId}")
                data.sent_message_len += len(recv_data)
                
        if mask & selectors.EVENT_WRITE:
            
            if data.content:
                if len(data.content) == 2:
                    print('sending the header')
                elif len(data.content) == 1:
                    print('sending the content')

                data.outb_reg = data.content.pop()
                #data.content = None

            if data.outb_reg:
                #print(f"Client : Sending {data.outb_reg}.")
                sent = sock.send(data.outb_reg)  # Should be ready to write
                data.outb_reg = data.outb_reg[sent:]

            elif not data.outb_reg :  # If not data to send, close the connection
                print(f"****** Closing connection {data.connId} ******")
                print("-------------------------------------------------------")
                self.selector.unregister(sock)
                sock.close()

    @staticmethod
    def prep_data (mode, outb_data, compression=False) :

        if mode == 'files' :  
            with open(outb_data, "rb") as f :
                proc_outb_data = f.read()

                header = str({
                    'mode' : mode,
                    'filename': outb_data}).encode()

            message_packet = types.SimpleNamespace(
                content= [header + b'$SPLIT$' + proc_outb_data ],
                outb_reg=b"",
                content_len=len(outb_data),
                sent_message_len=0,
                connId = rnd.randint(1,1000000)
                ) 
            

            # TODO: Add compression here
        elif mode in ['messages', 'coords', 'IR', 'robot_st', 'robot_end', 'dust'] :  

            proc_outb_data = outb_data.encode()

            header = str({
                    'mode' : mode,
                    'filename': None}).encode()

            message_packet = types.SimpleNamespace (
                    content = [header + b'$SPLIT$' + proc_outb_data ],
                    outb_reg=b"",
                    content_len=len(proc_outb_data),
                    sent_message_len=0,
                    connId = rnd.randint(1,1000000)
                    )   
        else:
            raise TypeError ('Mode must be either "file" or "message"')

        return message_packet


    def run (self):

        print("Client started running.")


        # Keeps polling to find new messages to send
        try:
            print('current message is', m)
        # while True:
            #TODO: By threading enable initiating new connections while the previous ones are running in the background
            
            # print("Waiting for new connections.")
            # print("-------------------------------------------------------")
            # print("Please enter your new connection in format shown below:")
            # inp_ = input ("Usage: <host(str)> <port(int)> <mode>\n")
            # print("----------")
            # message = input ("Please enter the message or the name of the file you want to send: \n")
            # print("----------")
            host= '127.0.0.1'
            #host, port, mode = inp_.split()[0], int(inp_.split()[1]), inp_.split()[2]
            proc_outb_data = self.prep_data (mode, m, compression=False)

            # initate a new connection
            self.initiate_connection(
                host=host,
                port=self.port,
                message_packet=proc_outb_data,
            )
            while self.selector.get_map(): # Keep sending till the connection is closed
                events = self.selector.select(timeout=1)
                if events:
                    for key, mask in events:
                        self.service_connection(key, mask)

        except KeyboardInterrupt:
            print("Caught keyboard interrupt, exiting")
        finally:
            self.selector.close()


def simulate (typ) :
        if typ == 'camera' :
            return 'files', [f'img/sim_camera/{num}.jpeg' for num in range(1,8)]
        elif typ == 'IR' :
            return 'IR', [str([[rnd.randint(-60,60) for _ in range(10)] for _ in range(100)])]
        elif typ == 'dust' :
            xs = list(range(20))
            es = np.random.normal(scale=6, size=len(xs))
            ys = np.sin(xs) + es
            return 'dust', [str(y) for y in ys]
            #return 'messages', [str(y) for y in list(range(20))]
        elif typ == 'coords':
            return 'coords', [str([rnd.randint(-60,60) for _ in range(7)])]
        elif typ == 'robot_st':
            return 'robot_st', ['robot commenced working']
        elif typ == 'robot_end':
            return 'robot_end', ['robot finished working']


if __name__ == '__main__' :
    port, simultor_mode = sys.argv[1], sys.argv[2]
    mode,  messages = simulate(simultor_mode)
    for m in messages :
        client = Client(port)
        client.run()


    