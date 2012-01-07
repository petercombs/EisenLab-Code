#===============================================================================
# Copyright (C) 2011 by Andrew Moffat
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#===============================================================================


import sys
import StringIO
import contextlib
import struct
import traceback
import inspect
import threading
import random
import socket
import json



import rlcompleter
try: import readline
except ImportError: readline = None


default_port = 1234







# this is used when we "exec" the code sent from the shell.  it essentially
# captures all stdout and returns it as a string.  we use a contextmanager
# to make absolutely sure that the real stdout is reassigned to sys.stdout
# no matter what happens 
@contextlib.contextmanager
def stdoutIO():
    stdout = StringIO.StringIO()
    old_out = sys.stdout
    sys.stdout = stdout
    yield stdout
    sys.stdout = old_out


class PortInUseException(Exception): pass




# actions that the shell can make to the server
COMMAND = "\x00"
AUTO_COMPLETE = "\x01"



def run_shell_server(f_globals, port):     
    auto_complete = rlcompleter.Completer(namespace=f_globals)
           
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    
    try: sock.bind(("localhost", port))
    except Exception, e:
        if isinstance(e, socket.error) and getattr(e, "errno", False) == 98:
            raise PortInUseException("%d in use" % port)
        else: raise
        
    sock.listen(100)   
    
    
    def run_repl(sock):
        
        def do_reply(data):
            data = json.dumps(data)
            out = struct.pack("i", len(data)) + data
            sock.send(out)

        while True:
            data = sock.recv(5)
            #print repr(data)
            
            if not data: return
            action, size = struct.unpack("!ci", data)

            data = sock.recv(size)
            if not data: return
            data = json.loads(data)
            
            
            if action == COMMAND:
                with stdoutIO() as stdout:
                    # lets us exec AND eval
                    try: exec compile(data, "<dummy>", "single") in f_globals, f_globals
                    except: print traceback.format_exc()
    
                out = stdout.getvalue()
                do_reply(out)
            
            elif action == AUTO_COMPLETE:
                ac = auto_complete.complete(*data)
                do_reply(ac)
            
    
    # our main loop for dispatching new connections
    while True:
        new_conn, addr = sock.accept()
        ct = threading.Thread(target=run_repl, args=(new_conn,))
        ct.daemon = True
        ct.start()
        






def open_shell(port):
    sock = socket.socket()
    sock.connect(("localhost", port))

    
    
    def do_request(action, data):
        data = json.dumps(data)
        data = action + struct.pack("!i", len(data)) + data
        sock.send(data)
    
        size = sock.recv(4)
        if not size: return
        size = struct.unpack("i", size)[0]
    
        if size:
            reply = sock.recv(size)
            return json.loads(reply)
        
            
        
    if readline:
        def auto_completer(*args):
            reply = do_request(AUTO_COMPLETE, args)
            return reply
        
        readline.set_completer(auto_completer)
        readline.parse_and_bind("tab: complete")
            
    
    prompt = "is:%d> " % port

    while True:
        try: line = raw_input(prompt)
        except EOFError: return
        
        reply = do_request(COMMAND, line)
        if reply is not None: print reply
            
            









# running from commandline, open a shell
if __name__ == "__main__":    
    try: port = int(sys.argv[1])
    except: port = default_port
    
    if readline:
        import os
        import atexit
        
        histfile = os.path.join(os.path.expanduser("~"), ".inspect_shell_history")
        try: readline.read_history_file(histfile)
        except IOError: pass
        atexit.register(readline.write_history_file, histfile)
    
    open_shell(port)
    
# it's being imported, run the shell server
else:
    caller = inspect.stack()[1][0]
    f_globals = caller.f_globals
    port = f_globals.get("inspect_shell_port", default_port)
    
    st = threading.Thread(target=run_shell_server, args=(f_globals, port))
    st.daemon = True
    st.start()
    