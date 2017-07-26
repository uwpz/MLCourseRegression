
df.ames = read_delim("./data/AmesHousing.txt", delim = "\t")
colnames(df.ames) = str_replace_all(colnames(df.ames), " ", "_")
summary(mutate_if(df.ames, is.character, as.factor))


hist(log(df.ames$Screen_Porch+1))



"C:/My/Stat/My"


ggplot(mpg, aes(displ, hwy)) +
  geom_point() +
  facet_wrap(~class)


summary(mpg)

library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)
ggplot(mpg, aes(cty, displ)) +
  geom_bin2d(bins = 100) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_fill_gradientn(colours=r, trans="log")

mtcars


l.split = split(1:nrow(mtcars), (1:nrow(mtcars)) %/% 10)
test = foreach(i = 1:length(l.split), .combine = c) %dopar% {
  tmp = mtcars[l.split[[i]],"mpg"]
}


df.tmp = mtcars %>% top_n(2) %>% select(carb)
top_n(df.tmp, 2)

count_(mtcars, "carb", sort=TRUE)$carb[1:2]

library(survival)
df.tmp = aml
fit <- survfit(Surv(time, status) ~ x, data = aml) 
plot(fit, lty = 2:3)


