######################################################
##### mosquito biting time from the Burkina Faso #####
######################################################

# data from https://malariajournal.biomedcentral.com/articles/10.1186/s12936-020-03538-5#Sec15

nia <- rbind(readxl::read_excel("data/Quantifying_mosquito_bites_BF.xlsx", sheet = "Mosquito data Niakore indoors")[,c("Village", "Location", "Hour", "Species", "Sex")],
             readxl::read_excel("data/Quantifying_mosquito_bites_BF.xlsx", sheet = "Mosquito data Niakore outdoors")[,c("Village", "Location", "Hour", "Species", "Sex")]) %>% subset(Sex == "F")

bt_data <- rbind(nia[,c("Village", "Location", "Hour", "Species")],
                 readxl::read_excel("data/Quantifying_mosquito_bites_BF.xlsx", sheet = "Mosquito data Toma indoors")[,c("Village", "Location", "Hour", "Species")],
                 readxl::read_excel("data/Quantifying_mosquito_bites_BF.xlsx", sheet = "Mosquito data Toma outdoors")[,c("Village", "Location", "Hour", "Species")]) %>% subset(Species == "An. gambiae sl")

bt_data$Hour <- factor(bt_data$Hour, levels = c("07pm-08pm", "08pm-09pm", "09pm-10pm", "10pm-11pm", "11pm-12am", "12am-01am",
                                                "01am-02am", "02am-03am", "03am-04am", "04am-05am", "05am-06am"))

bt_count_indoor <- bt_data %>% subset(Location == "IN")
bt_density_indoor <- bt_count_indoor %>% group_by(Hour) %>% summarise(count = n(),
                                                                      d = count/nrow(bt_count_indoor))

bt_count_outdoor <- bt_data %>% subset(Location == "OUT")
bt_density_outdoor <- bt_count_outdoor %>% group_by(Hour) %>% summarise(count = n(),
                                                                        d = count/nrow(bt_count_outdoor))

bt_density <- bt_data %>% group_by(Hour) %>% summarise(count = n(),
                                                       d = count/nrow(bt_data))

bt_density_outdoor$s_time <- c(seq(19,23), seq(0, 5))
bt_density_indoor$s_time <- c(seq(19,23), seq(0, 5))

bt_density_all <- rbind(bt_density_outdoor %>% mutate(Location = "Outdoor"),
                    bt_density_indoor %>% mutate(Location = "Indoor"))