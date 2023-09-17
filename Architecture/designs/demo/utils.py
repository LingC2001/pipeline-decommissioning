import datetime

def store_message (output, data) -> None:    
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

def write_file (output_dir, data) :
    file = open(output_dir, "ab")
    #conn.send("Filename received.".encode(FORMAT))
    """ Receiving the file data from the client. """
#    print(f"[RECV] Receiving the file data.")
    file.write(data)
#    conn.send("File data received".encode(FORMAT))
