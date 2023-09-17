import socket

# Client
SERVER_IP_ADDRESS = "127.0.0.1" # The IP address of the server
PORT_NO = 45678
MESSAGE = b'Hello Pipe! We are about to cut through you =))'

def client_send_message (
    protocol = 'TCP',
    receiver_ip = SERVER_IP_ADDRESS,
    port_no = PORT_NO,
    message = MESSAGE,
    ) -> None:

    """
    Sends a given message to the given server
    ...
    Arguments
    ----------
    protocol : str
        TCP/UDP
    receiver_ip : str
        Server's IP address
    port_no : int
        Server's listening port no.
    message : Unicode (encoded string)
        Message to be sent
    """

    if protocol == 'UDP':
        with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as clientSock:
            clientSock.connect((receiver_ip, port_no))
            clientSock.sendall(message)

    elif protocol == 'TCP':
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as clientSock:
            clientSock.connect((receiver_ip, port_no))
            clientSock.sendall(message)
    else:
        print (f'Invalid prtocol, choose between TCP and UDP')

    print('Data sent.')

if __name__ == '__main__' :
    client_send_message(
        protocol='TCP'
    )