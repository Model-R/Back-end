#le nomes
ary <- read.delim("./data/tax flora/nomes_ary.txt")
ary <- distinct(ary)###tinha muitas filas!
dim(ary)
ary <- mutate(ary,nome=paste(genero,epiteto))%>%arrange(nome.original)#cria coluna com nome
head(ary)

library(flora)
library(taxize)
library(taxonstand)

for(i in 1:length(ary$nome_spp)) ary$no.authors[i] <- remove.authors(ary$nome_spp[i])

flora.pack1 <- get.taxa(ary$nome, replace.synonyms = TRUE, suggest.names = TRUE,
                life.form = FALSE, habitat = FALSE, vernacular = FALSE,
                states = FALSE, establishment = FALSE, drop = c("infra.epiteth"),
                suggestion.distance = 0.9)

     ##write.table(flora.pack1,"florapack.txt")

     not.found2016 <- c("Algernonia pardina",
     "Annona amambayensis",
     "Archontophoenix cunninghamiana",
     "Capsicum parviflorum",
     "Cedrela balansae",
     "Celtis tala",
     "Chrysochlamys saldanhae",
     "Cnicothamnus lorentzii",
     "Cordiera stipulacea",
     "Daphnopsis americana",
     "Eugenia leitonii",
     "Gyrocarpus acuminatus",
     "Hieronyma oblonga",
     "Ilex congesta",
     "Inga tripa",
     "Jatropha palmatipartita",
     "Lonchocarpus nicou",
     "Luehea fiebrigii",
     "Luehea microcarpa",
     "Lycium cestroides",
     "Matayba leucodictya",
     "Miconia theaezans",
     "Morus alba",
     "Morus nigra",
     "Ocotea megaphylla",
     "Ormosia costulata",
     "Poincianella echinata",
     "Pouteria fragrans",
     "Roupala meisneri",
     "Ruprechtia crenata",
     "Sageretia elegans",
     "Schinus pearcei",
     "Sebastiania argutidens",
     "Sebastiania edwalliana",
     "Simira podocarpa",
     "Sloanea stipitata",
     "Tabebuia gemmiflora",
     "Talisia espiritosantensis",
     "Tovomita paniculata",
     "Trichilia trifolia",
     "Vitex laciniosa",
     "Ziziphus mistol")

     tpl.nf <-  TPL(not.found2016)
     tnrs.nf <- tnrs(not.found2016)
