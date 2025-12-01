# An attempt to create a GUI for the APS Gateways
import tkinter as tk
import threading
import numpy as np
import sys
import serial
import time
import re
import json
from collections import deque
from datetime import datetime
import wave
import pyaudio
import serial.tools.list_ports

def read_dict_from_json(filepath):
    try:
        with open(filepath, 'r') as f:
            data = json.load(f)
        return data
    except FileNotFoundError:
        print(f"Error: File not found at {filepath}")
        return None
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON format in {filepath}")
        return None

#TODO: State of Health?
# Logging (check?)
# Single Trial (schedule a node, make sure each node writes a file)
# Make regex for DIR{int}FIL{int}

# TODO: For 11/6/2025
# Schedule Trial
# Make Sure that Disconnect is working
# Make sure that initial node connect is working
# Move airmar thread to its own window

# node 1 \xF6\xDC\xBF\xBF\x19\x40
# node 2 \xEA\x3B\x0B\xE6\xC9\x64
baudrate = 115200

nodesmap = {}  # mapping of hyperion addr to nib ids
nodesmap = read_dict_from_json('nodesmap.json')

logfile = open("logfile.txt", "w")
logfile.write(f"M D Y: {datetime.now().strftime('%m %d %Y')}\n")
logfile.close()

nibreadlogfile = open("nibreadlogfile.txt", "w")
nibreadlogfile.write(f"M D Y: {datetime.now().strftime('%m %d %Y')}\n")
nibreadlogfile.close()

pps_to_wait = 3
time_to_play = 1

global_nodes = []

remain_open = True

class Gateway:
    def __init__(self, master, port, index):
        self.master = master
        self.index = index
        self.in_port = port[0]
        self.out_port = port[1]
        self.command_queue = deque()
        self.wait_period = 5
        self.play_period = 5
        print("in ser")
        self.in_ser = serial.Serial(self.in_port, baudrate, timeout=1)
        print("out ser")
        self.out_ser = serial.Serial(self.out_port, baudrate, timeout=1)

        print("in thread")
        self.in_thread = threading.Thread(target = self.in_thread)
        self.in_thread.start()
        print("out thread")
        self.out_thread = threading.Thread(target = self.out_thread)
        self.out_thread.start()

        # if (len(port) == 3):
        #     self.debug_port = port[2]
        #     self.debug_ser = serial.Serial(self.debug_port, baudrate, timeout=1)
        #     self.debug_thread = threading.Thread(target = self.debug_thread)
        #     self.debug_thread.start()

        self.nodes = []

        self.curnode = -1
        self.curhex = -1

    def parse_line(self, line):

        # print("Parsing: ", line)
        decode = line[3:-2].decode() # cut of I: and \r\n
        decode = decode
        if decode != "PPS DETECTED":
            print("Parsing: ", line)
            print("Decoded: ", decode)
        # print("Decoded: ", decode)
        # Possible messages:
        #startup
        #   I: USB: Baudrate detected: 115200\r\n
        #start scanning
        #   I: Start BLE Scanning\r\n
        #   I: BLE: Scanning Started!\r\n
        #   I: BLE: Connected: {6 byte adress} (random)\r\n
        #       I: BLE: Connected: EA:3B:0B:E6:C9:64 (random)\r\n
        #   I: BLE: Discovery Complete!\r\n
        #start ncpa mode
        #NIB status
        #schedule node to speak
        #ERROR
        #   I: UNKNOWN COMMAND\r\n
        if decode[0:4] == "BLE:":
            split = decode.split(": ")
            if split[1] == "Connected":
                hexstring = split[2][0:17]
                hs = hexstring.split(":")
                node = "".join(hs)
                self.nodes.append(node)
                global_nodes.append(node)
                print("Adding node: ", node)
                if node not in nodesmap:
                    self.master.create_radio_button(text=node)
                else: #TODO: check that a radio button isnt already in the set? i.e., its in nodesmap, but not already made
                    self.master.create_radio_button(text=node + ":" + nodesmap[node])
            elif split[1] == "Disconnected":
                hexstring = split[2][0:17]
                hs = hexstring.split(":")
                node = "".join(hs)
                self.nodes.remove(node)
                global_nodes.remove(node)
                print("Removing node: ", node)
                for i, rb in enumerate(self.master.radio_buttons):
                    if (not rb is None) and node in rb.cget("text"):
                        print(f"Found matching button: {node} and {rb.cget('text')}")
                        self.master.radio_buttons.pop(i)
                        rb.destroy()
                        break
            elif split[1] == "Discovery Complete!":
                print("discovery completed")
        elif decode[0:4] == "USB:":
            split = decode.split(": ")
            if split[1] == "INVALID HYP COMMAND!":
                print("Invalid hyp command!")
            elif baudrate != int(split[-1]):
                print(f"Baudrates dont match! {baudrate} vs {int(split[-1])}")
        elif "Sensor NIB read:" in decode:
            self.curnode = str("".join(decode[1:18].split(":")))
            print(f"Set curnode to {self.curnode}")
        elif decode == "UNKNOWN COMMAND":
            print("Error in sending code, unknown code")
        elif decode == "Start BLE Scanning":
            print("Start BLE Scanning")
        elif decode == "Write NIB":
            print("Write NIB")
        elif decode == "NIB Write ACK":
            print("NIB Write ACK")
        elif decode == "READ FROM NIB":
            print("READ FROM NIB")
        elif "CMD Write ACK" in decode:
            print("CMD Write ACK")
        elif decode == "Send CMD Flag":
            print("Send CMD Flag")
        elif decode == "Read SOH CMD":
            print("Read SOH CMD")
        elif decode == "Turn off PPS Mode":
            print("Turning off PPS Mode")
        elif decode == "Turn on PPS Mode":
            print("Turning on PPS Mode")
        elif decode == "Get Gateway SOH":
            print("Getting Gateway SOH")
        elif decode == "PPS DETECTED":
            print(end='')
        elif decode == "Error signing data":
            print("Error signing data")
        elif decode == "Write NIB - PPS":
            print("Write NIB - PPS")
        elif "PPS Timeout" in decode:
            print("PPS Timeout")
        elif ("Current" in decode) and ("Volt" in decode) and ("Mode" in decode):
            print("SOH power: ", decode)
        elif ("Sensor SOH read:" in decode):
            self.curhex = decode[1:18].replace(":", "")
        elif ("GPS" in decode) and ("Lat" in decode) and ("Lon" in decode):
            print("SOH gps: ", decode)
            status = decode.split(" ")[-1].split(":")[1]
            # print("Status is ", status)
            # print("Curhex is ", self.curhex)
            if self.curhex != -1:
                for rb in self.master.radio_buttons:
                    print(rb.cget("text"))
                    if (not rb is None) and (self.curhex in rb.cget("text")):
                        # print("Found matching button")
                        if status == 'V':
                            rb.configure(bg = "red")
                        elif status == 'A':
                            rb.configure(bg = "green")
            self.curhex = -1
        elif "File Count" in decode:
            print("SOH file: ", decode)
        elif self.curnode != -1 and self.curnode not in nodesmap: #TODO: Update nib_read response parsing
            try:
                line2 = str(bytearray.fromhex(line.decode()).decode())
                nibreadlogfile = open("nibreadlogfile.txt", "a")
                nibreadlogfile.write(self.curnode+":"+line2)
                nibreadlogfile.close()
                if ("DIR" in line2) and ("FIL" in line2):
                    # matches = re.findall("DIR\d*FIL\d*", line2)
                    matches = re.findall("\d+", line2)
                    if len(matches) == 2:
                        print(f"{self.curnode}: DIR {matches[0]} FIL {matches[1]}")
                        nib_read_return.append(matches[1])
                        #TODO: change later
                        if len(self.nodes) == len(nib_read_return):
                            print(f"All nodes reported back, latest fiels: {nib_read_return}")
                    else:
                        print("Wrong number of matches in Dir/Fil nib_read!")
                        return
                else:
                    # should look like "READY--NODEXXX!!!..." up to 32 bytes
                    matches = re.findall("\d{3}", line2)
                    if len(matches) == 1:
                        nodesmap[self.curnode] = matches[0]
                        print(f"Mapped {self.curnode} to {nodesmap[self.curnode]}")
                        print(nodesmap)
                        with open("nodesmap.json", "w", encoding="utf-8") as f: # save new nodesmap
                            json.dump(nodesmap, f, indent=4, ensure_ascii=False)
                        for rb in self.master.radio_buttons:
                            if (not rb is None) and rb.cget("text") == self.curnode:
                                rb.config(text = self.curnode + ":" + nodesmap[self.curnode])
                                break
                    else:
                        print("Wrong number of matches in ready nib_read!", line2)

                self.master.text += line2 + "\n"
                print("self.master.text: ", self.master.text)
                self.master.update_text()
                self.curnode = -1
            except Exception as e:
                print("Exception: ", e, ", with line: ", line)
                self.curnode = -1
        elif self.curnode != -1 and self.curnode in nodesmap:
            # a nib read after the fact!
            try:
                line2 = str(bytearray.fromhex(line.decode()).decode())
                nibreadlogfile = open("nibreadlogfile.txt", "a")
                nibreadlogfile.write(self.curnode + ":" + line2)
                nibreadlogfile.close()
                if ("DIR" in line2) and ("FIL" in line2):
                    # matches = re.findall("DIR\d*FIL\d*", line2)
                    matches = re.findall("\d+", line2)
                    if len(matches) == 2:
                        print(f"{self.curnode}: DIR {matches[0]} FIL {matches[1]}")
                        nib_read_return.append(matches[1])
                        #TODO: change later
                        if len(self.nodes) == len(nib_read_return):
                            print(f"All nodes reported back, latest fiels: {nib_read_return}")
                    else:
                        print("Wrong number of matches in Dir/Fil nib_read!")
                        return
                print("From nib read: ", line2)
                self.master.text += line2 + "\n"
                print("self.master.text: ", self.master.text)
                self.master.update_text()
                self.curnode = -1
            except Exception as e:
                print(f"Exception {e} in {line}")
                self.curnode = -1
        else:
            print("Unhandled line: ", line)

        return line

    def start_scanning(self, command):
        print("Sending command: ", command)
        # self.out_ser.write(command)
        self.command_queue.append(command)

    def start_ncpa_mode(self, command):
        for node in self.nodes:
            # cmd = command | bytes(f"\x00\x00\x00{node}\x00\x00\x00\x00", 'utf-8') # bitwise or to write in specific node adress
            #cmd = command | b"\x00\x00\x00"+node+b"\x00\x00\x00\x00"
            # cmd = bytes(map(lambda a,b: a | b, command, bytes(f"\x00\x00\x00{node}\x00\x00\x00\x00", 'utf-8')))
            cmd = command[0:3] + bytearray.fromhex(node) + command[-4:]
            print("Sending start_ncpa_mode command: ", cmd)
            # self.out_ser.write(cmd)
            self.command_queue.append(cmd)
            time.sleep(.5) # TODO: find better way to do this? After ACK, maybe?

    def stop_ncpa_mode(self, command):
        for node in self.nodes:
            cmd = command[0:3] + bytearray.fromhex(node) + command[-4:]
            print("Sending start_ncpa_mode command: ", cmd)
            # self.out_ser.write(cmd)
            self.command_queue.append(cmd)
            time.sleep(.5) # TODO: find better way to do this? After ACK, maybe?

    def read_soh(self, command):
        for node in self.nodes:
            cmd = command[0:3] + bytearray.fromhex(node) + command[-4:]
            # print("Bytearray.fromhex: ", bytearray.fromhex(node), bytearray.fromhex(node)[-6:])
            print("Sending read_soh command: ", cmd)
            # self.out_ser.write(cmd)
            self.command_queue.append(cmd)
            time.sleep(.5)

    def nib_read(self, command):
        for node in self.nodes:
            cmd = command[0:4] + bytearray.fromhex(node)
            print("Sending nib_read command: ", cmd)
            # self.out_ser.write(cmd)
            self.command_queue.append(cmd)
            time.sleep(1)  # TODO: find better way to do this? After ACK, maybe?


    def nib_write(self, node):
        nodeid = int(nodesmap[node])
        # nodeid = 100
        numseconds = pps_to_wait
        cmd = f"NIBWS{nodeid:04d}{numseconds:1d}00000000000000000000000000".encode()
        print("Sending nib_write command: ", cmd)
        # self.out_ser.write(cmd)
        self.command_queue.append(cmd)
        # time.sleep(numseconds+2) # + 1 is for wiggle room after
        # play_sound()

    def nib_write_silent(self):
        nodeid = 0
        numseconds = pps_to_wait
        cmd = f"NIBWS{nodeid:04d}{numseconds:1d}00000000000000000000000000".encode()
        print("Sending nib_write_silent command: ", cmd)
        # self.out_ser.write(cmd)
        self.command_queue.append(cmd)
        # time.sleep(numseconds+2) # + 1 is for wiggle room after

    def nib_write_correlate(self):
        nodeid = 0
        numseconds = pps_to_wait
        cmd = f"NIBWT{nodeid:04d}{numseconds:1d}00000000000000000000000000".encode()
        print("Sending nib_write_correlate command: ", cmd)
        # self.out_ser.write(cmd)
        self.command_queue.append(cmd)
        # time.sleep(numseconds+2) # + 1 is for wiggle room after

    def nib_write_loud(self):
        nodeid = 0
        numseconds = pps_to_wait
        cmd = f"NIBWA{nodeid:04d}{numseconds:1d}00000000000000000000000000".encode()
        print("Sending nib_write_loud command: ", cmd)
        # self.out_ser.write(cmd)
        self.command_queue.append(cmd)
        # time.sleep(numseconds+2) # + 1 is for wiggle room after

    def read_gateway_soh(self, command):
        print("Sending read_gateway_soh command: ", command)
        self.command_queue.append(command)
        time.sleep(.5)

    def start_pps(self, command):
        # for node in self.nodes:
        #     cmd = command[0:3] + bytearray.fromhex(node) + command[-4:]
        print("Sending start_pps command: ", command)
        self.command_queue.append(command)
        time.sleep(.5)

    def stop_pps(self, command):
        # for node in self.nodes:
        #     cmd = command[0:3] + bytearray.fromhex(node) + command[-4:]
        print("Sending stop_pps command: ", command)
        self.command_queue.append(command)
        time.sleep(.5)

    def in_thread(self):
        while remain_open:
            if(self.in_ser.inWaiting() > 0):
                line = self.in_ser.readline()
                logfile = open("logfile.txt", "a")
                logfile.write(f"{datetime.now().strftime('%H:%M:%S:%f')} Gateway {self.index}, reading {line.decode()[:-2]}\n")
                logfile.close()
                # print("Read line: ", line)
                parsed = self.parse_line(line)
            else:
                time.sleep(.1)

    def out_thread(self):
        while remain_open:
            if(len(self.command_queue) > 0):
                tmp_command = self.command_queue.popleft()
                logfile = open("logfile.txt", "a")
                logfile.write(f"{datetime.now().strftime('%H:%M:%S:%f')} Gateway {self.index}, sending {tmp_command}\n")
                logfile.close()
                self.out_ser.write(tmp_command)
            else:
                time.sleep(.1)

    def debug_thread(self):
        while remain_open:
            if(self.debug_ser.inWaiting() > 0):
                line = self.debug_ser.readline()
                logfile = open("logfile.txt", "a")
                logfile.write(f"{datetime.now().strftime('%H:%M:%S:%f')} Gateway {self.index}, debug {line.decode()[:-2]}\n")
                logfile.close()
                print(f"{datetime.now().strftime('%H:%M:%S:%f')} Gateway {self.index}, debug {line.decode()[:-2]}\n")
            else:
                time.sleep(.1)

nib_read_return = []



class ControlWindow(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        self.geometry("500x650")
        self.title("Control GUI")
        self.gateways = []
        self.center = tk.Frame(self, width=500, height=800)
        self.center.pack()
        self.frame1 = tk.Frame(self.center, width=250, height=800)
        self.frame1.pack(side="left")
        self.frame2 = tk.Frame(self.center, width=250, height=800)
        self.frame2.pack(side="left")
        self.create_button(self.frame1, "Start Scanning", self.start_scanning)
        self.create_button(self.frame1, "Read SOH", self.read_soh)
        self.create_button(self.frame1, "NIB Read", self.nib_read)
        self.create_button(self.frame1, "Gateway SOH", self.read_gateway_soh)
        self.create_button(self.frame2, "Start NCPA Mode", self.start_ncpa_mode)
        self.create_button(self.frame2, "Stop NCPA Mode", self.stop_ncpa_mode)
        self.create_button(self.frame2, "Start PPS", self.start_pps)
        self.create_button(self.frame2, "Stop PPS", self.stop_pps)
        self.schedule_trial_button = self.create_button(self.frame1, "Schedule Trial", self.schedule_trial)



    def create_button(self, frame, text="None", cmd = None):  # Create a button within the frame given
        button_width = 30
        button_height = 10

        self.button = tk.Button(frame, text=text, command=cmd)
        self.button.configure(height=button_height, width=button_width)
        self.button.pack()
        return self.button

    def add_gateways(self, ports):
        for i, port in enumerate(ports):
            self.gateways.append(GatewayWindow(master=self, port=port, index=i))
    def add_airmar(self):
        self.airmar = AirmarWindow(master=self)

    def start_scanning(self):
        scan_command = b'\x48\x59\x50\xFF\xFF\xFF\xFF\xFF\xFF\x09\xBB\xFF\xFF'
        for gateway in self.gateways:
            gateway.start_scanning(scan_command)
        return

    def start_ncpa_mode(self):
        blank_gps_only_mode_command = b'\x48\x59\x50\x00\x00\x00\x00\x00\x00\x03\x07\xFF\xFF'
        for gateway in self.gateways:
            gateway.start_ncpa_mode(blank_gps_only_mode_command)
        time.sleep(4)
        blank_start_ncpa_mode_command = b'\x48\x59\x50\x00\x00\x00\x00\x00\x00\x03\xAA\xFF\xFF' # missing \x00 bytes are written in later
        for gateway in self.gateways:
            gateway.start_ncpa_mode(blank_start_ncpa_mode_command)

        return

    def stop_ncpa_mode(self):
        blank_stop_ncpa_mode_command = b'\x48\x59\x50\x00\x00\x00\x00\x00\x00\x03\xBB\xFF\xFF'
        for gateway in self.gateways:
            gateway.stop_ncpa_mode(blank_stop_ncpa_mode_command)
        return

    def read_soh(self):
        blank_read_soh_command = b'\x48\x59\x50\x00\x00\x00\x00\x00\x00\x04\xFF\xFF\xFF'
        for gateway in self.gateways:
            gateway.read_soh(blank_read_soh_command)

    def nib_read(self):
        blank_nib_read_command = b'\x4E\x49\x42\x52\x00\x00\x00\x00\x00\x00'
        nib_read_return = []
        for gateway in self.gateways:
            gateway.nib_read(blank_nib_read_command)
        return

    def read_gateway_soh(self):
        blank_read_gateway_soh_command = b'\x48\x59\x50\xFF\xFF\xFF\xFF\xFF\xFF\x13\xFF\xFF\xFF'
        for gateway in self.gateways:
            gateway.read_gateway_soh(blank_read_gateway_soh_command)

    def start_pps(self):
        blank_start_pps_command = b'\x48\x59\x50\xFF\xFF\xFF\xFF\xFF\xFF\x11\xFF\xFF\xFF'
        for gateway in self.gateways:
            gateway.start_pps(blank_start_pps_command)

    def stop_pps(self):
        blank_start_pps_command = b'\x48\x59\x50\xFF\xFF\xFF\xFF\xFF\xFF\x12\xFF\xFF\xFF'
        for gateway in self.gateways:
            gateway.stop_pps(blank_start_pps_command)

    def schedule_trial(self):
        self.schedule_trial_button.configure(bg = 'green')
        self.update()
        for node in global_nodes:
            for gateway in self.gateways:
                gateway.select_node(node)
                gateway.schedule_node_direct(node)
            self.update()
            time.sleep(pps_to_wait + time_to_play + 9) # wait until its supposed to play, wait another second to play, add some wiggle room on top
            # INSERT HERE gateway.schedule_node_correlate()
            # INSERT HERE extra sleep for correlation response
            self.nib_read() # log responses
            time.sleep(1)
        print("Trial Done!")
        self.schedule_trial_button.configure(bg = 'white')

class AirmarWindow(tk.Toplevel):
    def __init__(self, master):
        tk.Toplevel.__init__(self, master)
        self.master = master
        self.geometry("500x500")
        self.frame = tk.Frame(self, width=500, height=500)
        self.frame.pack()
        self.title("Airmar Window")
        self.text = ""
        self.textArea = tk.Text(self.frame, height=500, width=500)
        self.textArea.pack()
        self.airmarcomport = "COM4"
        self.airmarbaudrate = 4800

        print("airmar thread")
        self.airmar_thread = threading.Thread(target=self.airmar_thread)
        self.airmar_thread.start()

    def format_dict(self, dict):
        outputstr = ""
        for key in dict:
            outputstr += key + str(dict[key]) + "\n"
        return outputstr

    def parse_nmea(self, line):
        split = line.split(",")
        type = split[0]
        checksum = split[-1]
        if type == "$HCHDT":
            #True heading, https://receiverhelp.trimble.com/alloy-gnss/en-us/nmea0183-messages-hdt.html
            if len(split) != 3:
                print("$HCHDT message is incorrect length")
                return -1
            heading = split[1]
            return heading
        elif type == "$YXXDR":
            #Transducer measurements
            return -1
        elif type == "$WIMWV":
            #WIMWV is standard wind speed and angle, http://www.nuovamarea.net/blog/wimwv
            if len(split) != 6:
                print("$WIMWV message is incorrect length")
                return -1
            windangle = split[1]
            ref = split[2]
            windspeed = split[3]
            windspeedunit = split[4]
            valid = split[5][0]
            return (windangle, ref, windspeed, windspeedunit, valid)
        elif type == "$GPZDA":
            #$GPZDA is UTC datetime, https://docs.novatel.com/OEM7/Content/Logs/GPZDA.htm
            if len(split) != 7:
                print("$GPZDA message is incorrect length")
                return -1
            utc = split[1]
            day = split[2]
            month = split[3]
            year = split[4]
            localzone = split[5] # ???
            localmin = split[6] # ??
            return (utc, day, month, year)
        elif type == "$WIMDA":
            #WIMDA is Meteorlogical Composite, http://www.nuovamarea.net/blog/wimda
            if len(split) != 21:
                print("$WIMDA message is incorrect length")
                return -1
            barompresmerc = split[1]
            barompresbar = split[3]
            airtempc = split[5]
            watertempc = split[7]
            relativehumidity = split[9]
            absolutehumidity = split[10]
            dewpointc = split[11]
            winddirtrue = split[13]
            winddirmagnetic = split[15]
            windspeedknots = split[17]
            windspeedmps = split[19]
            return (barompresmerc, barompresbar, airtempc, watertempc, relativehumidity, absolutehumidity, dewpointc, winddirtrue, winddirmagnetic, windspeedknots, windspeedmps)
        elif type == "$GPGGA":
            return -1
        elif type == "$GPVTG":
            return -1
        else:
            print("Uknown case: ", type)
            return -1

    def airmar_thread(self):
        self.airmarcomport = serial.Serial(self.airmarcomport, self.airmarbaudrate, timeout=1)
        self.msgdict = {}
        while remain_open:
            if (self.airmarcomport.inWaiting() > 0):
                block = self.airmarcomport.read(4096)
                self.airmarlogfile = open("airmarlogfile.txt", "ab")
                self.airmarlogfile.write(block)
                self.airmarlogfile.close()
                blockconv = block.decode()
                for line in blockconv.split("\n"):
                    # print(line)
                    if len(line) > 0 and line[0] =='$':
                        msgtype = line.split(",")[0]
                        self.msgdict[msgtype] = line.split(",")[1:]
                self.textArea.delete("1.0", tk.END)
                self.textArea.insert(tk.END, self.format_dict(self.msgdict))
            else:
                time.sleep(.1)

class GatewayWindow(tk.Toplevel):
    def __init__(self, master, port, index):
        tk.Toplevel.__init__(self, master)
        self.index = index
        self.geometry("500x800")
        self.frame = tk.Frame(self, width=500, height=800)
        self.frame.pack()
        self.radio_frame = tk.Frame(self.frame, width=5100, height=300)
        self.radio_frame.pack()
        self.radio_buttons = []
        self.title(f"Gateway GUI {port[0]}, {port[1]}")
        self.port = port
        self.gateway = Gateway(self, self.port, self.index)
        self.text = " "
        self.textarea = tk.Text(self, height=500, width=300)
        self.textarea.pack()

        self.create_button("Schedule Node", self.schedule_node_wrapper)
        self.create_button("Silent Schedule", self.schedule_node_silent_wrapper)
        self.create_button("Loud Schedule", self.schedule_node_loud_wrapper)
        self.create_button("Correlate", self.schedule_node_correlate_wrapper)

    def create_radio_button(self, text="None"):
        rb = tk.Radiobutton(self.radio_frame, text=text, variable = stringvar, value=text)
        rb.pack(side=tk.TOP, ipady=5)
        stringvar.set("")
        self.radio_buttons.append(rb)

    def create_button(self, text="None", cmd = None):  # Create a button within the frame given
        button_width = 30
        button_height = 8

        self.button = tk.Button(self.frame, text=text, command=cmd)
        self.button.configure(height=button_height, width=button_width)
        self.button.pack()

    def update_text(self):
        self.textarea.delete("1.0", tk.END)
        self.textarea.insert(tk.END, self.text)

    def start_scanning(self, command):
        self.gateway.start_scanning(command)
        return

    def start_ncpa_mode(self, command):
        self.gateway.start_ncpa_mode(command)

    def stop_ncpa_mode(self, command):
        self.gateway.stop_ncpa_mode(command)

    def nib_read(self, command):
        self.text = ""
        self.gateway.nib_read(command)

    def read_soh(self, command):
        self.gateway.read_soh(command)

    def schedule_node_wrapper(self):
        for gateway in self.master.gateways:
            gateway.schedule_node()
        time.sleep(pps_to_wait + 2)  # + 1 is for wiggle room after
    def schedule_node_silent_wrapper(self):
        for gateway in self.master.gateways:
            gateway.schedule_node_silent()
        time.sleep(pps_to_wait + 2)  # + 1 is for wiggle room after
    def schedule_node_correlate_wrapper(self):
        for gateway in self.master.gateways:
            gateway.schedule_node_correlate()
        time.sleep(pps_to_wait + 2)  # + 1 is for wiggle room after
    def schedule_node_loud_wrapper(self):
        for gateway in self.master.gateways:
            gateway.schedule_node_all()
        time.sleep(pps_to_wait + 2)  # + 1 is for wiggle room after

    def schedule_node(self):
        if stringvar.get() != '':
            if ":" in stringvar.get():
                nodeid = stringvar.get().split(":")[0]
                print("Scheduling node: ", nodeid)
                self.gateway.nib_write(nodeid)
            else:
                print(f": not in stringvar: {stringvar.get()}")
        else:
            print("Select node to schedule!")

    def select_node(self, nodeid):
        for rb in self.radio_buttons:
            if (not rb is None) and (nodeid in rb.cget("text")):
                stringvar.set(rb.cget("text"))

    def schedule_node_direct(self, nodeid):
        self.gateway.nib_write(nodeid)

    def schedule_node_silent(self):
        self.gateway.nib_write_silent()

    def schedule_node_correlate(self):
        self.gateway.nib_write_correlate()

    def schedule_node_all(self):
        self.gateway.nib_write_loud()

    def read_gateway_soh(self, command):
        self.gateway.read_gateway_soh(command)

    def start_pps(self, command):
        self.gateway.start_pps(command)

    def stop_pps(self, command):
        self.gateway.stop_pps(command)

def make_chirp():
    chirp = wave.Wave_write("./chirp.wav")
    numseconds = 3
    chirp.setparams((1, 2, 48000, 48000*numseconds, 'NONE', 'not compressed'))
    # numchans, sampwidth, samplingrate, numsamps, comptype, compname = chirp.getparams()
    numchans, sampwidth, samplingrate, numsamps, comptype, compname = (1, 2, 48000, 48000*numseconds, 'NONE', 'not compressed')
    t = np.linspace(0, numseconds, numsamps, endpoint=False)
    audio_data = ((5000 * (np.sin(2 * np.pi * 440 * t) + 0.5 * np.sin(2 * np.pi * 720 * t) + 0.33 * np.sin(
        2 * np.pi * 1050 * t)))).astype(np.int16)
    # print(audio_data)
    # audio_data = np.concatenate((audio_data, audio_data, audio_data, audio_data, audio_data))
    # numsamps = numsamps*5
    chirp.writeframes(audio_data.tobytes())
    chirp.close()
    onechan = audio_data

def play_sound():
    f = wave.open("./chirp.wav")
    p = pyaudio.PyAudio()
    chunk = 1024
    stream = p.open(format=p.get_format_from_width(f.getsampwidth()),
                    channels=f.getnchannels(),
                    rate=f.getframerate(),
                    output=True)
    # read data
    data = f.readframes(chunk)

    # play stream
    while data:
        stream.write(data)
        data = f.readframes(chunk)

        # stop stream
    stream.stop_stream()
    stream.close()

    # close PyAudio
    p.terminate()

def on_closing():
    global remain_open
    remain_open = False

make_chirp()

app = ControlWindow()
app.protocol("WM_DELETE_WINDOW", on_closing)
stringvar = tk.StringVar(app, "")
# Port should be the read from nrf and the port for MDASN (output)

# comports = serial.tools.list_ports.comports()
# for port, desc, hwid in sorted(comports):
#         print("{}: {} [{}]".format(port, desc, hwid))

# the read for nrf is the higher one, the lower one is debug?
# ports = [("COM8", "COM10")] # A list of tuples of (input, output) for each of the desired gateways
# ports = [("COM11", "COM10", "COM12")]
ports = [("COM23", "COM25", "COM24"), ("COM11", "COM10", "COM12")] # ("gateway1input", "gateway1output", "gateway1debug"), ("gateway2input"...) ON HARLEYS COMPUTER
#ports = [("COM10", "COM11", ""), ("COM3", "COM6", "")] # ON RUGGED
# ports = [("/dev/ttyACM1", "/dev/ttyACM3")]
app.add_gateways(ports)
# app.add_airmar()

app.mainloop()