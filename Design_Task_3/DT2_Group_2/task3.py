
import random
import numpy as np
import math

class Queue:
    def __init__(self):
        self.list = []

    def push(self,item):
        self.list.insert(0,item)

    def pop(self):
        return self.list.pop()

    def isEmpty(self):
        return len(self.list) == 0
    

def get_succ(rotations:list):
    succ = list()
    
    copy_of = rotations[:]
    copy_of.append(0)
    succ.append(copy_of)
    
    copy_of = rotations[:]
    copy_of.append(-1)
    succ.append(copy_of)
    
    copy_of = rotations[:]
    copy_of.append(1)
    succ.append(copy_of)

    return succ


def is_goal(rotations:list, angle_deg, use_init=False):

    z0 = angle_deg

    if (use_init):
        d_init = 1

        if (angle_deg < 0):
            d_init = - 1

        z0 = angle_deg - (d_init * 90)

    z = z0

    for i, d_i in enumerate(rotations):
        if (d_i != 0):
            z = z - (d_i * np.rad2deg(np.arctan(2 ** (-i))))

    if (abs(z) <= 0.2):
        return True, abs(z)
    return False, abs(z)
        



def find_rotation_strings(angle_deg, use_init=False):

    frontier = Queue()
    initial_state = list()
    visited_states = list()
    frontier.push(initial_state)
    
    valid_rotations = list()
    rotations_errors = list()

    iteration_limit = 100000
    i = 0

    

    while(i < iteration_limit):

        current_state = frontier.pop()

        if (len(valid_rotations) != 0 and len(current_state) > len(valid_rotations[0])):
            return valid_rotations, rotations_errors

        goal, error = is_goal(current_state, angle_deg, use_init)

        if (goal):
            valid_rotations.append(current_state)
            rotations_errors.append(error)


        
        if (current_state not in visited_states):
            visited_states.append(current_state)
            succ = get_succ(current_state)
            for next_rotations in succ:
                if (next_rotations not in visited_states):
                    frontier.push(next_rotations)
        i += 1
    
    return False
    



def main():
    angle_deg = -33.75
    use_init = False

    print("Finding rotation strings for angle: {}°\n".format(angle_deg))

    rotations, errors = find_rotation_strings(angle_deg, use_init)

    if not rotations:
        print("Could't find solution")
        return

    print("Shortest rotation strings and their errors are: ")
    best_zero_count = 0
    best_error = errors[0]
    best_string = rotations[0] 

    for i, rotation_string in enumerate(rotations):
        print("{}, error = {}".format(rotation_string, errors[i]))
        zero_count = rotation_string.count(0)
        error = errors[i]
        if (zero_count >= best_zero_count):
            if (zero_count == best_zero_count and error < best_error):
                best_zero_count = zero_count
                best_error = error
                best_string = rotation_string


    print("\nBest string (least rotations / highest zero count and lowest error) is : {}".format(best_string))


    cordic_gain = 1
    for i in range(len(best_string)):
        if (best_string[i] != 0):
            cordic_gain *= np.sqrt(1 + 2 ** (-2 * i))

    print("CORDIC gain for the best string: {}".format(cordic_gain))

    x = random.uniform(-20, 20)
    y = random.uniform(-20, 20)

    print("--------------------------------------------------------------------------------------------------------------------------------------------------------------")
    print("Testing the rotation with values (x = {}, y = {})\n".format(x, y))

    if (use_init):
        d_init = 1
        if (angle_deg < 0):
            d_init = -1

        x0 = -d_init * y
        y0 = d_init * x

        x_i = x0
        y_i = y0
    else:
        x_i = x
        y_i = y

    for i in range(len(best_string)):
        x_n = x_i - (best_string[i] * y_i * (2 ** (-i)))
        y_n = y_i + (best_string[i] * x_i * (2 ** (-i)))

        x_i = x_n
        y_i = y_n

    x_prime = x_i / cordic_gain
    y_prime = y_i / cordic_gain
    print("Result of the rotation using CORDIC with best string:{} (x_cord = {},{} y_cord = {}) ".format(" "*3, x_prime, " "*4, y_prime))

    # Reference using Givens
    x_prime_ref = x * np.cos(np.deg2rad(angle_deg)) - y * np.sin(np.deg2rad(angle_deg))
    y_prime_ref = y * np.cos(np.deg2rad(angle_deg)) + x * np.sin(np.deg2rad(angle_deg))

    print("Result of the rotation using Givens:{} (x_giv = {},{} y_giv = {}) \n".format(" "*20,x_prime_ref, " "*5, y_prime_ref))
   

    u = np.array([x_prime, y_prime])
    v = np.array([x_prime_ref, y_prime_ref])
    dot_product = np.dot(u, v)
    magnitude_u = np.linalg.norm(u)
    magnitude_v = np.linalg.norm(v)
    angle_radians = math.acos(dot_product / (magnitude_u * magnitude_v))
    angle_degrees = math.degrees(angle_radians)

    print("Angle between CORDIC rotated and Givens rotated vector: {}°".format(angle_degrees))
    print("--------------------------------------------------------------------------------------------------------------------------------------------------------------")





if __name__ == "__main__":
    main() 