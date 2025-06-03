import os
import shutil
import sys
from itertools import combinations
import multiprocessing
import re
import time



def worker_function(lock, MEND_file,process_id):

    temperatures=[1.5,2]
    #temperature = [x for x in temperatures]
    temperature = [x - 2.01 for x in temperatures]
    with lock:
        print(f"Process {process_id} is modifying the file.")
        time.sleep(30)
        # Define the replacement value
        import re
        shutil.copyfile(MEND_file, 'MEND_namelist.nml')

        # Create a list to store the modified lines
        modified_lines = []
        filename='MEND_namelist.nml'
        with open(filename, 'r') as file:
            for line in file:
                if 'Dir_Output  =' in line:
                    folder_prefix=line.split('userio/')[1].strip().replace("'", '')
                if '    sSTP_delta  = 0.0' in line:
                    modified_line = '    sSTP_delta  ='+str(temperature[process_id])
                    modified_lines.append(modified_line)
                        
                else:
                    modified_lines.append(line)
                # Modify the lines as needed
            modified_lines[42] = "    Dir_Output  = 'userio/"+folder_prefix+str(process_id)+"'\n"
            folder = modified_lines[42].split('=')[1].strip().replace("'", '').replace('/', '\\')
        # Write the modified lines back to the file
        with open(filename, 'w') as file:
            file.writelines(modified_lines)  
            print(f"Process {process_id} finished modifying the file.")
        with open('userio/inp/tair_rh_rad_hourly_ambient.in', 'r') as infile:
            lines = infile.readlines()
        header = lines[0]
        data_lines = lines[1:]
        # Process the lines
        new_lines = [header.strip()]
        for line in data_lines:
            values = line.strip().split(",")
            # Assuming we are adding 10 to the second column
            values[0] = str(float(values[0]) + temperature[process_id])
            new_lines.append(','.join(values))

        # Write the new file
        with open('userio/inp/tair_rh_rad_hourly.in', 'w') as outfile:
            for line in new_lines:
                outfile.write(line + '\n')

    if folder:
        os.makedirs(folder, exist_ok=True)

    shutil.copy(filename, os.path.join(folder, 'MEND_namelist.nml'))
    os.system("python pipeline.py test")
    

if __name__ == "__main__":
    MEND_file='MEND_namelist_normal.nml'

    num_processes = 2

    lock = multiprocessing.Lock()


    processes = []
    for i in range(num_processes):
        process = multiprocessing.Process(target=worker_function, args=(lock, MEND_file, i))
        processes.append(process)
        try:
            process.start()
        except KeyboardInterrupt:
            print("Terminating main process.")
            process.terminate()
            process.join()
            sys.exit(1)

    for process in processes:
        process.join()




