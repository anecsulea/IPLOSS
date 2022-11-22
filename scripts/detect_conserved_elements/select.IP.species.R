###############################################################

pathResults="../../results/whole_genome_alignments/"

###############################################################

library(ape)

###############################################################

tree=read.tree(file=paste(pathResults, "avian_366/Birds366_tree.txt", sep=""))
species=tree$tip.label

###############################################################

withIP=c("Chauna_torquata", "Anser_cygnoid", "Cairina_moschata", "Asarcornis_scutulata", "Anas_zonorhyncha", "Anas_platyrhynchos", "Anas_platyrhynchos_platyrhynchos", "Anseranas_semipalmata", "Pauxi_pauxi", "Penelope_pileata", "Casuarius_casuarius", "Dromaius_novaehollandiae", "Eudromia_elegans", "Nothocercus_nigrocapillus", "Nothocercus_julius", "Nothoprocta_perdicaria", "Nothoprocta_ornata", "Nothoprocta_pentlandii", "Rhea_americana", "Rhea_pennata", "Struthio_camelus", "Tinamus_guttatus")

###############################################################

print(paste(withIP, collapse=","))

###############################################################

tree.withIP=keep.tip(tree, withIP)

write.tree(tree.withIP, paste(pathResults, "species_with_IP/tree.txt",sep=""))

###############################################################
