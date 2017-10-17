library (RColorBrewer)
library (gplots) 

dat1 <- read.table('/Users/jrteubl/Desktop/NEW/col1a1_matrix.txt', header = TRUE, sep = '\t', row.names= 1)
dat_mat1 <- data.matrix(dat1)

mammals <- c('Little.Brown.Bat', 'Large.Fruit.Bat', 'Nine.banded.Armadillo', 'Tasmanian.Devil', 'Gray.Short.tailed.Opossum',
      'Tammar.Wallaby', 'European.Hedgehog', 'Shrew','American.Pika', 'European.Rabbit',
      'Platypus', 'African.Bush.Elephant', 'Lesser.Hedgehog.Tenrec', 'Northern.Treeshrew','Bottlenose.Dolphin')
birds <- c('Mallard', 'Chicken', 'Turkey', 'Collared.Flycatcher', 'Zebra.Finch')
fishes <- c('Japanese.Rice.Fish', 'Southern.Platyfish', 'Atlantic.Cod', 'Three.spined.Stickleback', 'Spotted.Gar',
          'Nile.Tilapia', 'Japanese.Puffer', 'Green.Spotted.Puffer')
reptiles <- c('Carolina.Anole', 'Chinese.Soft.shell.Turtle')
amphibians <- c('Western.Clawed.Frog')
artidactyls <- c('Cattle', 'Sheep', 'Wild.Boar', 'Alpaca', 'Horse')
carnivores <- c('Giant.Panda', 'Dog', 'Cat', 'European.Polecat')
primates <- c('Marmoset', 'Green.Monkey', 'Western.Gorilla', 'Human', 'Rhesus.Macaque', 'Northern.White.cheeked.Gibbon',
              'Common.Chimpanzee', 'Olive.Baboon', 'Philippine.Tarsier')
rodents <- c('Guinea.Pig', 'Ords.Kangaroo.Rat', 'Striped.Gopher', 'Mouse', 'Rat')

crb <- list()

for (i in colnames(dat1)){
  if (i %in% mammals){
    crb <- append(crb, 'Thistle')
  } else if (i %in% birds){
    crb <- append(crb, 'LightGreen')
  } else if (i %in% fishes){
    crb <- append(crb, 'Gold')
  } else if (i %in% reptiles){
    crb <- append(crb, 'Turquoise')
  } else if (i %in% amphibians){
    crb <- append(crb, 'MediumPurple')
  } else if (i %in% artidactyls){
    crb <- append(crb, 'DarkSlateBlue')
  } else if (i %in% carnivores){
    crb <- append(crb, 'LightCoral')
  } else if (i %in% primates){
    crb <- append(crb, 'SkyBlue')
  } else if (i %in% rodents){
    crb <- append(crb, 'SeaGreen')
  } else sprintf("Could not find %s", i)
}
crb <- unlist(crb)

legendNames <- c("mammal(other)", "bird", "fish", "reptile", "amphibian", "ungulate", "carnivore", "primates", "rodents")
colorNames <- c('Thistle', 'LightGreen', 'Gold', 'Turquoise', 'MediumPurple', 'DarkSlateBlue', 'LightCoral', 'SkyBlue', 'SeaGreen')
pdf(file = '/Users/jrteubl/Desktop/NEW/heatmap_col1a1.pdf')
par(cex.main=1)
par(xpd = T, mar = par()$mar + c(0,3,0,7))
heatmap.2(dat_mat1, dendrogram = 'row', col= brewer.pal(9,"BuGn"),
          scale = "none",trace = "none", symm = TRUE, offsetRow = -0.5, offsetCol = 0,
          margins=c(14,18), cexRow=.7, cexCol = .7,
          keysize = 1, key.title = NA, key.xlab = NA, key.ylab = NA,
          main = "Ratio of shared peptides\nin COL1A1", RowSideColors = crb)
legend("topleft",
       legend = legendNames, # category labels
       col = colorNames,  # color key
       lty= 1,             # line style
       lwd = 5,            # line width
       cex = .5,
       bty = "n",
       ncol=2
       
       
)

dev.off()

