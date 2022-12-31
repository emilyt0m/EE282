library(ggplot2)

#Sequence length distribution
length_short <- read.table("length_short.txt", header = FALSE)
length_long <- read.table("length_long.txt", header = FALSE)
colnames(length_short) <- c('name', 'length', 'assembly' )
colnames(length_long) <- c('name', 'length', 'assembly' )
length_short$log_length <- log(length_short$length)
length_long$log_length <- log(length_long$length)

plot_length_short <- ggplot(data = length_short, aes (x=log_length)) + geom_histogram(bins=200) + labs(title = 'Log Sequence Length, <100kb')
plot_length_short

plot_length_long <- ggplot(data = length_long, aes (x=log_length)) + geom_histogram(bins=7) + labs(title = 'Log Sequence Length, >=100kb')
plot_length_long

#GC% distribution
gc_short <- read.table("gc_short.txt", header = FALSE)
gc_long <- read.table("gc_long.txt", header = FALSE)
colnames(gc_short) <- c('name', 'GC', 'assembly' )
colnames(gc_long) <- c('name', 'GC', 'assembly' )

plot_gc_short <- ggplot(data = gc_short, aes (GC)) + geom_histogram(bins=200) + labs(title = 'GC %, <100kb')
plot_gc_short

plot_gc_long <- ggplot(data = gc_long, aes (x=GC)) + geom_histogram(bins=20) + labs(title = 'GC %, >=100kb')
plot_gc_long
