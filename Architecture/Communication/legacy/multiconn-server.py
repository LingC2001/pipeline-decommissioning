#!/usr/bin/env python3

import sys
import socket
import selectors
import types
from utils import *

sel = selectors.DefaultSelector()


def accept_wrapper(sock):
    conn, addr = sock.accept()  # Should be ready to read
    print(f"Accepted connection from {addr}")
    conn.setblocking(False)
    data = types.SimpleNamespace(
        addr=addr,
        inb=b"",
        outb=b""
        )
    events = selectors.EVENT_READ | selectors.EVENT_WRITE
    sel.register(conn, events, data=data)


def service_connection(key, mask):
    sock = key.fileobj
    data = key.data
    if mask & selectors.EVENT_READ:
        recv_data = sock.recv(1024)  # Should be ready to read
        print(f"Server received : {recv_data}")

        if recv_data:

            if mode == 'file' :
                write_file (filename, recv_data)
            elif mode == 'message' :
                store_message ('messages.txt', recv_data)
            else:
                raise TypeError ('The 4th keword arg. "mode", must be file or message.')

            # data.outb = recv_data
            # # data.outb += (b'recieved' + recv_data)
            # data.inb += recv_data

        if not recv_data : # empty message found in the buffer
            print(f"Closing connection to {data.addr}")
            sel.unregister(sock)
            sock.close()

    if mask & selectors.EVENT_WRITE:
        if data.outb:
            print(f"Echoing {data.outb!r} to {data.addr}")
            sent = sock.send(data.outb)  # Should be ready to write
            data.outb = data.outb[sent:]
            # sock.sendall(data.outb)
            


if len(sys.argv) != 5:
    print(f"Usage: {sys.argv[0]} <host> <port> <mode> <filename/message>")
    sys.exit(1)

# TODO: Get file name from client itself.

host, port, mode, filename = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]

# Initiate a TCP socket
lsock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
lsock.bind((host, port))
lsock.listen()
print(f"Listening on {(host, port)}")
lsock.setblocking(False)
sel.register(lsock, selectors.EVENT_READ, data=None)

# TODO : timeout conenction

try:
    while True:
        events = sel.select(timeout=None)
        for key, mask in events:
            if key.data is None: 
                accept_wrapper(key.fileobj) # Start a new connection
            else:
                service_connection(key, mask)
except KeyboardInterrupt:
    print("Caught keyboard interrupt, exiting")
finally:
    sel.close()
