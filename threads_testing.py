import threading
import time

def worker_task(name, delay):
    print(f"Task {name} started...")
    time.sleep(delay)
    print(f"Task {name} finished!")

# Create the threads
thread1 = threading.Thread(target=worker_task, args=("A", 2))
thread2 = threading.Thread(target=worker_task, args=("B", 3))

# Start the threads
thread1.start()
thread2.start()

# Wait for both threads to finish before moving the main script forward
thread1.join()
thread2.join()

print("All tasks complete!")
