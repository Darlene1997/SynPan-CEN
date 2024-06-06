library(RColorBrewer)
library(ggplot2)
library(Cairo)
args<-commandArgs(T)  #NIP_SL177  15
prefix = args[1]
variable <- args[2]

monomer_match <- paste(prefix, paste0("12chr_for_plot_ED", variable, ".txt"), sep="_")
colorrange <- paste(prefix, paste0("12chr_ColorRange_ED", variable, ".txt"), sep="_")

######### monomer_match
df_match <- read.table(monomer_match)
names(df_match) <- c("chr", "id", "start", "end", "y", "y_pair", "material", "pair", "family")
df_match$family <- factor(df_match$family, levels=c("F1", "SF1","F2", "SF2","F3", "SF3","F4", "SF4","F5", "SF5","F6", "SF6","F7", "SF7","F8", "SF8","F9", "SF9","F10", "SF10","F11", "SF11","F12", "SF12","F13", "SF13","F14", "SF14","F15", "SF15", "LTR", "gene", "gap", "blank", "plus", "minus", "SZ-22", "RIRE", "otherLTR", "DNAtrans", "box", "Interval"))

######### pair color range
df_colorrange <- read.table(colorrange)
names(df_colorrange ) <- c("chr", "x", "y", "pair", "order", "ED")
df_colorrange$ED <- factor(df_colorrange$ED, levels=c("0","1", "2","3", "4","5", "6","7", "8","9", "10","11", "12","13", "14","15", "LTR"))

plot_all <- ggplot() +
  geom_rect(data=df_match[df_match$family != "box", ], mapping=aes(xmin=start, xmax=end, ymin=y-0.3, ymax=y+0.3, fill=family))+
  geom_rect(data=df_match[df_match$family == "box", ], mapping=aes(xmin=start, xmax=end, ymin=y-0.3, ymax=y+0.3), fill="transparent", color="black", linewidth=0.05)+
  geom_polygon(data=df_colorrange, mapping=aes(x=x , y=y , group=pair , fill=ED))+
  scale_fill_manual(values = c("F1"="#3f66a1", "SF1"="#3f66a1","F2"="#9D5427", "SF2"="#9D5427","F3"="#85A7CC", "SF3"="#85A7CC","F4"="#D0AF62", "SF4"="#D0AF62","F5"="#D67C1B", "SF5"="#D67C1B","F6"="#3897c5", "SF6"="#3897c5","F7"="#a874b5", "SF7"="#a874b5","F8"="#7ec0b4", "SF8"="#7ec0b4","F9"="#6DAB30", "SF9"="#6DAB30","F10"="#096858", "SF10"="#096858","F11"="#d66c54", "SF11"="#d66c54","F12"="#b13f73", "SF12"="#b13f73","F13"="#afa7d8", "SF13"="#afa7d8","F14"="#6566a9", "SF14"="#6566a9","F15"="#893f8b", "SF15"="#893f8b", "LTR"="#A9A9A9", "gene"="black", "gap"="red", "blank"="#D3D3D3", "0"="#004FAF","1"="#006995","2"="#00837B","3"="#009E60","4"="#00B846","5"="#00D32B","6"="#3DFF00","7"="#57FF00","8"="#72FF00","9"="#8CFF00", "10"= "#A7FF00", "11"= "#C1FF00", "12"= "#DBFF00", "13"= "#FFD300", "14"= "#FFB800", "15"= "#FF9E00", "plus"= "red", "minus"= "blue", "SZ-22"= "#5B86A8", "RIRE"= "#AC7DB3", "otherLTR"= "#F4B617", "DNAtrans"= "#C1D6E7", "box"= "transparent", "Interval"= "#262626"))+
  scale_color_manual(values = c("F1"="#3f66a1", "SF1"="#3f66a1","F2"="#9D5427", "SF2"="#9D5427","F3"="#85A7CC", "SF3"="#85A7CC","F4"="#D0AF62", "SF4"="#D0AF62","F5"="#D67C1B", "SF5"="#D67C1B","F6"="#3897c5", "SF6"="#3897c5","F7"="#a874b5", "SF7"="#a874b5","F8"="#7ec0b4", "SF8"="#7ec0b4","F9"="#6DAB30", "SF9"="#6DAB30","F10"="#096858", "SF10"="#096858","F11"="#d66c54", "SF11"="#d66c54","F12"="#b13f73", "SF12"="#b13f73","F13"="#afa7d8", "SF13"="#afa7d8","F14"="#6566a9", "SF14"="#6566a9","F15"="#893f8b", "SF15"="#893f8b", "LTR"="#A9A9A9", "gene"="black", "gap"="red", "blank"="#D3D3D3", "0"="#004FAF","1"="#006995","2"="#00837B","3"="#009E60","4"="#00B846","5"="#00D32B","6"="#3DFF00","7"="#57FF00","8"="#72FF00","9"="#8CFF00", "10"= "#A7FF00", "11"= "#C1FF00", "12"= "#DBFF00","13"= "#FFD300", "14"= "#FFB800", "15"= "#FF9E00", "plus"= "red", "minus"= "blue", "SZ-22"= "#5B86A8", "RIRE"= "#AC7DB3", "otherLTR"= "#F4B617", "DNAtrans"= "#C1D6E7", "box"= "black", "Interval"= "#262626"))+
  facet_wrap(~chr, ncol = 1)+ scale_x_continuous(limits = c(0, 1500000), breaks = seq(0, 1500000, 100000))+
  theme(strip.background = element_blank(), panel.background = element_blank(), legend.position = "none")+
  theme(strip.text = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_text(size=10))+
  ggtitle(prefix)

output <- paste(prefix, paste0("matchList_12chr_pair_ED", variable, ".pdf"), sep="_")

CairoPDF(file=output, width=10, height = 8)
plot_all
dev.off()

