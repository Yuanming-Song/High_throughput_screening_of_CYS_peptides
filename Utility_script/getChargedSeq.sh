#!/bin/zsh

# Define common amino acid one-letter codes
amino_acids=(A R N D C Q E G H I L K M F P S T W Y V)

# Function to calculate the charge of a sequence
calculate_charge() {
    local sequence=$1
    local charge=0
    for ((i=1; i<=${#sequence}; i++)); do
        amino=$(echo "$sequence" | cut -c$i)
        case $amino in
            R|K|H) charge=$((charge + 1)) ;;  # Positively charged amino acids
            D|E) charge=$((charge - 1)) ;;    # Negatively charged amino acids
        esac
    done
    echo $charge
}

# Generate sequences and check charges
max_length=2  # Define the maximum peptide length
for length in {2..$max_length}; do
    echo "Generating sequences of length $length..."
    for position in {1..$length}; do
        echo "Position $position"
        # Generate sequences with C at the specified position
        for seq_comb in ${(C)amino_acids}; do
            sequence=" "
            for ((i=1; i<=length; i++)); do
                if ((i == position)); then
                    sequence+="C"
                else
                    sequence+="${seq_comb[$i]}"
                fi
            done

            # Calculate the charge
            charge=$(calculate_charge "$sequence")

            # Print sequence if charge is not zero
            if ((charge != 0)); then
                echo "$sequence"
            fi
        done
    done
done