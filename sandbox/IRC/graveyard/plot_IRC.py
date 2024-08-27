import matplotlib.pyplot as plt

# Function to read the .log file and extract IRC data including counting directions accurately
def extract_irc_data_and_count_directions(log_file_path):
    steps = []
    energies = []
    forward_count = 0
    reverse_count = 0
    current_direction = None

    with open(log_file_path, 'r') as file:
        for line in file:
            if "Point Number  1 in FORWARD path direction" in line:
                if current_direction != "FORWARD":  # To avoid recounting if multiple lines mention the direction
                    current_direction = "FORWARD"
            elif "Point Number  1 in REVERSE path direction" in line:
                if current_direction != "REVERSE":  # Switch direction only if it changes
                    current_direction = "REVERSE"
            elif "Summary of reaction path following" in line:
                break  # Stop reading further as we've reached the summary section

            if current_direction == "FORWARD" and "Point Number" in line:
                forward_count += 1
            elif current_direction == "REVERSE" and "Point Number" in line:
                reverse_count += 1

    # Since the counts include the 'Point Number: 1' line for both directions, we don't add 1 for the starting point here
    total_points = forward_count + reverse_count + 1
    
    return forward_count, reverse_count, total_points

# Function to plot the IRC energy profile remains unchanged
# This function will not be called in this snippet as it's focused on counting points

def plot_irc(steps, energies):
    plt.figure(figsize=(10, 6))
    plt.plot(steps, energies, marker='o', linestyle='-', color='blue')
    plt.title('IRC Energy Profile')
    plt.xlabel('Step Number')
    plt.ylabel('Energy (a.u.)')
    plt.grid(True)
    plt.show()


# Example usage:
log_file_path = 'H-N-IRC.log'  # Update this to the path of your .log file
forward_count, reverse_count, total_points = extract_irc_data_and_count_directions(log_file_path)

print(f"Forward points: {forward_count}")
print(f"Reverse points: {reverse_count}")
print(f"Total number of points: {total_points}")

plot_irc(steps, energies)


