import bugv #BuGVTrack BuGV

input_dir = "/public/home/jiajb/src/lib/python/bugv/"

bed_file = input_dir + "example/test.bed"
bam_file = input_dir + "test.RNA.bam"
sRNA_bam_file = input_dir + "test.sRNA.bam"

#gff_file = "/public/home/jiajb/lib_data/genome/core/ath/Araport11/a.gtf"
#gfr = bugff.GffRead(gff_file)
#df = gfr.load_region_data(["chr1", 3500, 6000, "+"])

bgv = bugv.BuGV()
bgv.add_track("RulerTrack",  config = {"track_height": 1, "ruler.ruler_length": 500})
bgv.add_track("BedTrack", file_name = bed_file, config = {"track_height": 1})
bgv.add_track("RNASashimiplot", file_name=bam_file, config={"track_height": 4, 
                                                            "ylim_style": "change_less",
                                                            "plot_junc": True,
                                                            "read_specific_strand": True,
                                                            #"ylim_up_down_ratio": 0.5,
                                                            "rnasashimiplot": {"line_width_based_on_count": True}
                                                           })
bgv.add_track("sRNAPointTrack", file_name=sRNA_bam_file, config={"track_height": 4,
                                                                "read_obj_kwargs": {"fileformat": "sim_seq_count"}})
bgv.add_track("SpaceTrack", config={"track_height": 0})
bgv.add_track("XaxisTrack", config={"track_height": 1})

bgv.load_track_config({"y_zero_line.color": "black"})
bgv.plot(region=["chr1", 2, 10000, "+"], select_features=None, dpi=600)
bgv.fig.savefig("/public/home/jiajb//test.png")