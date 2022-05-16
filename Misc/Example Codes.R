
# creating datasets ------------------------

circles <- read.table("Misc/circle_dataset.txt", head = F, sep = ",")
colnames(circles) <- c("x", "y")

moons <- read.table("Misc/moons_dataset.txt", head = F, sep = ",")
colnames(circles) <- c("x", "y")

left_moon <- read.table("Misc/left_moon.txt", head = T, sep = " ")[,1:2]
right_moon <- read.table("Misc/right_moon.txt", head = T, sep = " ")[,1:2]

variables <- c(rep("Left", nrow(left_moon)), rep("Right", nrow(left_moon)))

moons_qual_quant <- cbind(rbind(left_moon, right_moon), variables)
colnames(moons_qual_quant) <- c("x", "y", "category")
moons_qual_quant$category <- as.factor(moons_qual_quant$category)

single_moon <- read.table("Misc/single_moon_dataset.txt", head = F, sep = ",")
colnames(single_moon) <- c("x", "y")

usethis::use_data(single_moon, overwrite = TRUE)
