# Communication Protocols:

### What protocols can we use to establish and maintain communication between our server & client?

## 1. UDP

- User Datagram Protocol is connection-less protocol which is suitable for applications that require efficient communication that doesn't have to worry about packet loss.

> *Connectionless protocols* send signal without checking whether a receiver is available or available to receive the signal.

### **Comparison**:

1. Speed advantage over TCP
2. Lower Reliability

## 2. TCP

- A connection-oriented protocol, it relies on a server in a passive open state. A passive open server listens for any client trying to connect with it.

> The *connection* is established via a three-way handshake. The client sends a synchronization request, the server sends back an acknowledgment, and the client returns a synchronization acknowledgment in response. .

### **Comparison**:

1. Slower than UDP
2. Performs connection and packet checks

## 1. Multi-connections:

- For multi-connections we use selector library, which allows us to do I/O multiplexing. This allows for I/O completion checks on multiple sockets and hence our system will not be blocked by any socket (I/O operation). 

> Alternatively other concurrency methods such as async I/O could be used.