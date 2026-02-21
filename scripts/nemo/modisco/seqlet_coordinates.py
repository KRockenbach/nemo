import os, sys, h5py
import pandas as pd

# df associating motif names with TF family names
motif_family_df = pd.read_csv("../data/motifs/motif_families.tsv", delimiter="\t", index_col=None, header=0)

pattern_groups = ['pos_patterns', 'neg_patterns']


def create_df(organism):

    data_dict = {"seq": [], "family": [], "start": [], "end": [], "q_val": []}

    modisco_dir = f"../results/nemo/{organism}/masked_graphpart/modisco"

    for seq in ["promoter", "terminator"]:
        for segment_start in range(-500, 500, 100):

            segment_end = segment_start + 100
            segment_descriptor = f"{seq}_%{segment_start}to{segment_end}%"
            modisco_h5py = os.path.join(modisco_dir, f"{segment_descriptor}_modisco_results.h5")
            with h5py.File(modisco_h5py, 'r') as modisco_results:

                for name in pattern_groups:
                    if name not in modisco_results.keys():
                        continue
                    metacluster = modisco_results[name]
                    key = lambda x: int(x[0].split("_")[-1])

                    for pattern_name, pattern in sorted(metacluster.items(), key=key):
                        starts = pattern['seqlets']['start'][:][:]
                        ends = pattern['seqlets']['end'][:][:]
                        get_pos = lambda x: segment_start + x
                        starts = list(map(get_pos, starts))
                        ends = list(map(get_pos, ends))
                        # load tomtom_df for current pattern
                        tomtom_file = os.path.join(modisco_dir, segment_descriptor,"tomtom", f"{name}.{pattern_name}.tomtom.tsv")
                        tomtom_df = pd.read_csv(tomtom_file, delimiter="\t", index_col=None, header=0, skip_blank_lines=True, comment="#")

                        for motif in tomtom_df.loc[:,"Target_ID"].tolist():
                            if motif in motif_family_df.loc[:,"MOTIF"].tolist():
                                family = motif_family_df.loc[motif_family_df.loc[:,"MOTIF"]==motif,"FAMILY"].values[0]
                                q_val = tomtom_df.loc[tomtom_df.loc[:,"Target_ID"] == motif,"q-value"].values[0]

                                data_dict["start"].extend(starts)
                                data_dict["end"].extend(ends)
                                data_dict["seq"].extend([seq for _ in range(len(starts))])
                                data_dict["family"].extend([family for _ in range(len(starts))])
                                data_dict["q_val"].extend([q_val for _ in range(len(starts))])

    out_path = os.path.join(modisco_dir, "seqlet_coordinates.tsv")
    out_df = pd.DataFrame(data_dict)
    out_df.to_csv(out_path, sep='\t', header=True, index=False)


if __name__ == "__main__":
    for organism in ["Bnapus", "Athaliana"]:
        create_df(organism)


