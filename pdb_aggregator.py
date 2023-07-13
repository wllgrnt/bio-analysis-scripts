# for all pdb files in a directory, calculate the average error for each chain
# across residues.
import glob
import os


def calculate_average_error(file_name: str, chain="C") -> float:
    with open(file_name) as pdb_file:
        pdb = pdb_file.readlines()[1:-3]

    residues = [-1]
    errors = [-1.0]

    for atom in pdb:
        if atom[21] == chain:
            if int(atom[23:26]) != residues[-1]:
                residues.append(int(atom[23:26]))
                errors.append(float(atom[61:66]))
            else:
                pass

    errors = errors[1:]
    average_score = sum(errors) / len(errors)
    return average_score


if __name__ == "__main__":
    input_path = "230713"
    input_folder_paths = f"input_folder/{input_path}/f*"
    output_path = f"output_folder/{input_path}_scores.csv"
    input_folders = glob.glob(input_folder_paths)
    scores = []
    for input_folder in input_folders:
        pdb_files = os.path.join(input_folder, "f*.pdb")
        for pdb_file in glob.glob(pdb_files):
            print(f"calculating score for {pdb_file}")
            score = calculate_average_error(pdb_file)
            scores.append((pdb_file, score))

    # now we have the score for each file, just gotta make a csv.
    with open(output_path, "w") as output_file:
        output_file.write("pdb_file,score\n")
        for score in scores:
            output_file.write(f"{score[0]},{score[1]}\n")

    # expecting 5 scores and 45 folders, so 225 scores
    assert len(scores) == 225
