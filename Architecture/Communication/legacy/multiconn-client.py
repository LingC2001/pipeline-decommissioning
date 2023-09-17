import sys
import socket
import selectors
import types

sel = selectors.DefaultSelector()


def start_connections(host, port, num_conns):
    server_addr = (host, port)
    for i in range(0, num_conns):
        connid = i + 1
        print(f"Starting connection {connid} to {server_addr}")
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.setblocking(False)
        sock.connect(server_addr)
        #sock.connect_ex(server_addr)
        events = selectors.EVENT_READ | selectors.EVENT_WRITE
        data = types.SimpleNamespace(
            connid=connid,
            #data_len=sum(len(m) for m in out_data),
            recv_total=0,
            send_register=out_data,
            outb=b"",
        )
        sel.register(sock, events, data=data)


def service_connection(key, mask):

    sock = key.fileobj
    data = key.data

    if mask & selectors.EVENT_READ:
        #print(f"send register : {data.send_register}")
        recv_data = sock.recv(1024)  # Should be ready to read
        if recv_data:
            print(f"Client : Received {recv_data!r} from connection {data.connid}")
            data.recv_total += len(recv_data)
        
        #if not recv_data or data.recv_total == data.msg_total: # stop connection if empty message is received or total received is the same as total sent
        # print(f"Closing connection {data.connid}")
        # sel.unregister(sock)
        # sock.close()
            

    if mask & selectors.EVENT_WRITE:

        if not data.outb and data.send_register:
            data.outb = data.send_register
            data.send_register = None

        if data.outb:
            print(f"Client : Sending {data.outb!r} to connection {data.connid}")
            sent = sock.send(data.outb)  # Should be ready to write
            data.outb = data.outb[sent:]
            # sock.sendall(data.outb)
        
        elif not data.send_register :  # If not data to send, close the connection
            print(f"Closing connection {data.connid}")
            sel.unregister(sock)
            sock.close()


#TODO: Handle argument types here

if len(sys.argv) != 6:
    print(f"Usage: {sys.argv[0]} <host> <port> <num_connections> <mode> <file / message (list of comma separated strings)>")
    sys.exit(1)

host, port, num_conns = sys.argv[1:4]

if sys.argv[4] == 'file' :  
    with open(sys.argv[5], "rb") as f :
        out_data = f.read()
elif sys.argv[4] == 'message' :  
    out_data = [mess.encode() for mess in sys.argv[5].split(",")]
else:
    raise TypeError ('The 4th keword arg. "mode", must be file or message.')


start_connections(host, int(port), int(num_conns))

try:
    while True:
        events = sel.select(timeout=1)
        if events:
            for key, mask in events:
                service_connection(key, mask)
        # Check for a socket being monitored to continue.
        if not sel.get_map():
            break
        
except KeyboardInterrupt:
    print("Caught keyboard interrupt, exiting")
finally:
    sel.close()


