import socket
import datetime

SERVER_IP_ADDRESS = "127.0.0.1" # Accept traffic from all over internet
PORT_NO = 45678
SIGNALS_DB = '../Databases/commands.txt'

def server_st_listening(
    protocol = 'TCP',
    receiver_ip = SERVER_IP_ADDRESS,
    port_no = PORT_NO,
    store_message = True
    ) -> None:
    """
    Sends a given message to the given server
    ----------
    Arguments
    ----------
    protocol : str
        TCP/UDP
    receiver_ip : str
        Listening IP address(es)
    port_no : int
        Server's listening port no.
    """

    if protocol == 'UDP':
        with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as serverSock:
            # Binding the server to IP/port to listen to 
            serverSock.bind ((receiver_ip, port_no))
        
            try :
                while True:
                    data, addr = serverSock.recvfrom(1024) # blocks execution
                    
                    if store_message:
                        keep_message(data.decode())
                    print(f'data: {data}')
                    # acknowledgement
                    serverSock.sendto(b'recieved.', addr)

            except KeyboardInterrupt :
                print(f"Server : Connection Closed.")

    elif protocol == 'TCP':
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as serverSock:
            serverSock.bind((receiver_ip, port_no))
            serverSock.listen()
            conn, addr = serverSock.accept() # blocks execution till a connection is received
            # TODO: enable multiple connections
            with conn:
                while True:
                    data = conn.recv(1024)
                    print(data)
                    if not data:
                        break
                    else:
                        if store_message:
                            keep_message(data.decode()) 
                    # acknowledgement
                    conn.sendall(b'recieved')

        print(f"Server : Connection Closed.")

    else:
        print (f'Invalid prtocol, choose between TCP and UDP')



def keep_message (data, output=SIGNALS_DB) -> None:    
    """
    Writes the given message in "output".
    ----------
    Arguments
    ----------
    data : str
        Data to be stored.
    output : str 
        Output directory
    """
    with open(output, 'a') as f:
        f.writelines(f'{datetime.datetime.now()}, {data}\n')
    print('Data Stored')



if __name__ == '__main__' :
    server_st_listening(
        protocol='TCP'
    )