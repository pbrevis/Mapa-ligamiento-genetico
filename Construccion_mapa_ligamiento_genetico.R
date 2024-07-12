##Cargando paquetes de R
library(polymapR)
library(hexbin)
library(dplyr)

rm(list = ls()) # Limpiando el ambiente global

# Guardando la ubicación del directorio de trabajo,
# una vez que fuese definida manualmente
directorio <- getwd()
head(directorio)

##Leyendo archivo *.csv con datos simulados de genotipo
##en población tetraploide segregante bi-parental F1
##con x=5 grupos de ligamiento

dosage_df <- read.csv("tetraploid_dosage.csv",
                        stringsAsFactors = FALSE,
                        row.names = 1)
                        #primera columna tiene nombres de marcadores

#Vista rápida al df
glimpse(dosage_df)

#Confirmando que los datos fueron ingresados a un 'data frame'
class(dosage_df)

#polymapR requiere los datos en forma de matriz
tetraploid_pop <- as.matrix(dosage_df)

#Confirmando que los datos ahora están en forma de matriz
class(tetraploid_pop)

#Vista rápida a las primeras columnas (1-a-6) de la matriz de datos
head(tetraploid_pop[,1:6])


#Dimensiones de la matriz: filas (marcadores) x columnas (individuos)
#m=3000 marcadores
#n=207 individuos F1 + 2 padres
dim(tetraploid_pop)

#Chequeando que proporciones de segregación corresponden a lo esperado
#según dosis de los padres y herencia tetrasómica
f1checked <- checkF1(dosage_matrix = tetraploid_pop,
                     parent1 = "P1", parent2 = "P2",
                     F1 = colnames(tetraploid_pop)[3:ncol(tetraploid_pop)],
                     polysomic = TRUE, disomic = FALSE, mixed = FALSE,
                     ploidy = 4)
                     #Especie se considera autotetraploide


#Revisando resultados de desviación de la distribución de probabilidad.
#Altos valores (~1) de 'qall_mult' y 'qall_weights' indican consistencia
#entre dosis de los padres y segregación de marcadores en la población F1
glimpse(f1checked$checked_F1)

#check class
class(f1checked$checked_F1)

#check dimensions
dim(f1checked$checked_F1)

#Varios pasos a continuación para filtrar marcadores, dejando afuera
#aquellos que se desvían de proporciones de segregación esperadas
all_markers <- f1checked$checked_F1

#Resumen estadístico de 'qall_mult' y 'qall_weights'
summary(all_markers$qall_weights)
summary(all_markers$qall_mult)
quantile(all_markers$qall_weights, probs = c(.01, 0.1), na.rm = TRUE)

nonskewed_markers <- filter(all_markers, qall_weights > 0.85)
                     # sólo conservar marcadores cuyo valor qall_weights > 0.85

#Resumen estadístico de 'qall_mult' y 'qall_weights' dentro del subgrupo de
#marcadores seleccionados
summary(nonskewed_markers$qall_weights)
summary(nonskewed_markers$qall_mult)

nonskewed_marker_names <- nonskewed_markers$MarkerName

#Vista previa: nombre de marcadores confiables
head(nonskewed_marker_names)

#Número de marcadores confiables
length(nonskewed_marker_names)

#Porcentaje de marcadores seleccionados
100*(length(nonskewed_marker_names) / length(all_markers$MarkerName))

#Nuevo dataframe, excluyendo marcadores no confiables (skewed)
#Conservando 2731 de 3000 marcadores (~91%)
dosage_df2 <- dosage_df %>% filter(row.names(dosage_df) %in% nonskewed_marker_names)

#polymapR requiere los datos en forma de matriz
#nueva matriz tetraploid_pop2 ahora sólo con marcadores confiables
tetraploid_pop2 <- as.matrix(dosage_df2)

#Confirmando que los datos ahora están en forma de matriz
class(tetraploid_pop2)

#PCA de la dosis de marcadores, para detectar posibles individuos problemáticos
PCA_progeny(dosage_matrix = tetraploid_pop2, 
            highlight = list(c("P1", "P2")), 
            colors = "red")

#Opcional: eliminar individuo extremo por su nombre
#tetraploid_pop2 <- subset(tetraploid_pop2, select=-c(X16401_N019))


#Creando resumen sobre calidad de datos
markerdatasum <- marker_data_summary(dosage_matrix = tetraploid_pop2,
                           ploidy = 4,
                           pairing = "random",
                           parent1 = "P1",
                           parent2 = "P2",
                           progeny_incompat_cutoff = 0.05)

#Histograma de dosis en los padres (todos los marcadores, n=2731)
pq_before_convert <- parental_quantities(dosage_matrix = tetraploid_pop2, 
                                         las = 2)

#Convirtiendo tipos de segregación a su dosis más simple, y eliminando marcadores
#no segregantes para facilitar computación
segregating_data <- convert_marker_dosages(dosage_matrix = tetraploid_pop2, ploidy = 4)

#Histograma de dosis en los padres, luego de convertir tipos de segregación
pq_after_convert <- parental_quantities(dosage_matrix = segregating_data)

#Control de calidad de datos: eliminando marcadores (filas) con problemas
screened_data <- screen_for_NA_values(dosage_matrix = segregating_data, 
                                      margin = 1, # opción 1 indica marcadores
                                      cutoff =  0.1, # tolerancia de 10% NA
                                      print.removed = FALSE) 

#Control de calidad de datos: eliminando individuos (columnas) con problemas
screened_data2 <- screen_for_NA_values(dosage_matrix = screened_data, 
                                       margin = 2, # opción 2 indica individuos                                       
                                       cutoff = 0.1, # tolerancia de 10% NA

                                       print.removed = FALSE)

#Revisando individuos duplicados
screened_data3 <- screen_for_duplicate_individuals(dosage_matrix = screened_data2, 
                                                   cutoff = 0.95, 
                                                   plot_cor = TRUE)

#Revisando marcadores duplicados
screened_data4 <- screen_for_duplicate_markers(dosage_matrix = screened_data3)

#Opcional: identificar marcadores confiables en base al tamaño del bin
#Formula M/(4NLy)
#reliable.markers <- names(which(sapply(screened_data4$bin_list,length) >= 6))
#reliable_data <- screened_data4$filtered_dosage_matrix[reliable.markers,]

#Extrayendo la matriz de dosis, elemento al interior de 'screened_data4', que
#es producto de la función 'screen_for_duplicate_markers'
filtered_data <- screened_data4$filtered_dosage_matrix

#Histograma de dosis en los padres, luego de eliminar marcadores duplicados
#Esta gráfica muestra el número final de marcadores a usar en el desarrollo
#del mapa genético
pq_screened_data <- parental_quantities(dosage_matrix = filtered_data)



##INICIO DEL MAPEO

##Definiendo grupos de ligamiento (cromosomas) en base a marcadores simplex x nulliplex
##(tipo 1x0) del primer padre (P1)
##En especie estudiada, número de cromosomas = 7; homólogos = 28
SN_SN_P1 <- linkage(dosage_matrix = filtered_data, 
                    markertype1 = c(1,0), # marcador tipo 1x0
                    parent1 = "P1",
                    parent2 = "P2",
                    which_parent = 1, #seed parent
                    ploidy = 4,
                    pairing = "random")
#Vista rápida
glimpse(SN_SN_P1)

SN_SN_P1_linkage <- SN_SN_P1
SN_SN_P1_linkage$marker_pair <- paste(SN_SN_P1_linkage$marker_a,
                                     SN_SN_P1_linkage$marker_b)
SN_SN_P1_linkage %>% group_by(phase) %>% summarise(count = n_distinct(marker_pair))

#Gráfica de frecuencia de recombinación (r) vs. valores LOD
r_LOD_plot(linkage_df = SN_SN_P1, r_max = 0.5)

#Chequeando desviaciones en la relación esperada entre r y LOD
#en los marcadores tipo 1x0 del Padre 1
P1deviations <- SNSN_LOD_deviations(linkage_df = SN_SN_P1,
                                    ploidy = 4,
                                    N = ncol(filtered_data) - 2, #no. individuos F1
                                    alpha = c(0.05,0.2),
                                    plot_expected = TRUE,
                                    phase="coupling")


#Resumen
summary(P1deviations)

#Paso opcional: eliminar valores extremos que se desviían de la relación esperada
#entre r y LOD
#SN_SN_P1b <- SN_SN_P1[SN_SN_P1$phase == "coupling",][-which(P1deviations > 0.2),]

#Gráfica de frecuencia de recombinación (r) vs. valores LOD
#r_LOD_plot(linkage_df = SN_SN_P1b, r_max = 0.5)

#Agrupando marcadores tipo 1x0 del Padre 1 para identificar homólogos
P1_homologues <- cluster_SN_markers(linkage_df = SN_SN_P1, 
                                    LOD_sequence = seq(3, 9, 0.25),
                                    #rango LOD de 3 a 9, en pasos de 0.25
                                    LG_number = 5, #número de cromosomas
                                    ploidy = 4, #tetraploide
                                    parentname = "P1", #seed parent
                                    plot_network = FALSE,
                                    plot_clust_size = FALSE)

#Chequeando manualmente número de clusters, según valor LOD=3.5
P1_hom_LOD3.5 <- P1_homologues[["3.5"]] #LOD = 3.5
t3.5 <- table(P1_hom_LOD3.5$cluster)
print(paste("Número de clusters:",length(t3.5)))

#Número de marcadores por cluster ante un LOD = 3.5
t3.5[order(as.numeric(names(t3.5)))]

#Chequeando manualmente número de clusters, según valor LOD=5
P1_hom_LOD5 <- P1_homologues[["5"]] #LOD = 5
t5 <- table(P1_hom_LOD5$cluster)
print(paste("Número de clusters:",length(t5)))

#Número de marcadores por cluster ante un LOD = 5
t5[order(as.numeric(names(t5)))]

#Definiendo estructura de los grupos de ligamiento, en base a
#datos de marcadores tipo 1x0 del Padre 1
LGHomDf_P1 <- define_LG_structure(cluster_list = P1_homologues, 
                                  LOD_chm = 3.5, # LODs seleccionados
                                  LOD_hom = 5,   # en paso previo
                                  LG_number = 5)

#Vista previa a marcadores tipo 1x0 mapeados
head(LGHomDf_P1)

#Resumen de marcadores tipo 1x0 mapeados
summary(LGHomDf_P1)

#Vista a los resultados numéricos, función 'define_LG_structure' en Padre 1
table(LGHomDf_P1$LG, LGHomDf_P1$homologue)



## Opcional: Seleccionar sólo pares de marcadores en acoplamiento
# SN_SN_P1_coupl <- SN_SN_P1[SN_SN_P1$phase == "coupling",] 
# 
# #Agrupando los marcadores en fase de acoplamiento, usando valores altos de LOD
# P1_homologues_1 <- cluster_SN_markers(linkage_df = SN_SN_P1_coupl, 
#                                       LOD_sequence = c(3:12), 
#                                       LG_number = 41,
#                                       ploidy = 4,
#                                       parentname = "P1",
#                                       plot_network = FALSE,
#                                       plot_clust_size = FALSE)



#Calculando ligamiento entre marcadores tipo simplex (1x0) y duplex (2x0)
#del Padre 1
SN_DN_P1 <- linkage(dosage_matrix = filtered_data, 
                    markertype1 = c(1,0),
                    markertype2 = c(2,0),
                    which_parent = 1,
                    ploidy = 4,
                    pairing = "random")

#Usando nuevas asociaciones para establecer "puentes" entre grupos de
#ligamiento previamente identificados
LGHomDf_P1_1 <- bridgeHomologues(cluster_stack = P1_homologues[["6"]], 
                                 linkage_df = SN_DN_P1, 
                                 LOD_threshold = 4, 
                                 automatic_clustering = TRUE, 
                                 LG_number = 41,
                                 parentname = "P1")

#Vista a los resultados numéricos, función 'bridgeHomologues' en Padre 1
table(LGHomDf_P1_1$LG, LGHomDf_P1_1$homologue)






##Repitiendo todos los pasos de mapeo, ahora con Padre 2 (pollen parent)

##Ligamiento entre marcadores simplex (tipo 1x0)
SN_SN_P2 <- linkage(dosage_matrix = filtered_data, 
                    markertype1 = c(1,0),
                    parent1 = "P1",
                    parent2 = "P2",
                    which_parent = 2,
                    ploidy = 4,
                    pairing = "random")

#Vista rápida
glimpse(SN_SN_P2)

SN_SN_P2_linkage <- SN_SN_P2
SN_SN_P2_linkage$marker_pair <- paste(SN_SN_P2_linkage$marker_a,
                                      SN_SN_P2_linkage$marker_b)
SN_SN_P2_linkage %>% group_by(phase) %>% summarise(count = n_distinct(marker_pair))

#Gráfica de frecuencia de recombinación (r) vs. valores LOD
r_LOD_plot(linkage_df = SN_SN_P2, r_max = 0.5)

SN_SN_P2_coupl <- SN_SN_P2[SN_SN_P2$phase == "coupling",] 
                  # Sólo usando marcadores en fase de acoplamiento, m=6995

#Agrupando marcadores tipo 1x0 del Padre 2 para identificar homólogos
P2_homologues <- cluster_SN_markers(linkage_df = SN_SN_P2_coupl, 
                                    LOD_sequence = c(3:12), 
                                    LG_number = 5,
                                    ploidy = 4,
                                    parentname = "P2",
                                    plot_network = FALSE,
                                    plot_clust_size = FALSE)

#Chequeando manualmente número de clusters, según valor LOD=4
P2_hom_LOD4 <- P2_homologues[["4"]] #LOD = 4
t4 <- table(P2_hom_LOD4$cluster)
print(paste("Número de clusters:",length(t4)))

#Número de marcadores por cluster ante un LOD = 4
t4[order(as.numeric(names(t4)))]

#Calculando ligamiento entre marcadores tipo simplex (1x0) y duplex (2x0)
#del Padre 2
SN_DN_P2 <- linkage(dosage_matrix = filtered_data, 
                    markertype1 = c(1,0),
                    markertype2 = c(2,0),
                    which_parent = 2,
                    ploidy = 4,
                    pairing = "random")

LGHomDf_P2 <- bridgeHomologues(cluster_stack = P2_homologues[["6"]], 
                               linkage_df = SN_DN_P2, 
                               LOD_threshold = 4, 
                               automatic_clustering = TRUE, 
                               LG_number = 5,
                               parentname = "P2")

table(LGHomDf_P2$LG,LGHomDf_P2$homologue)

# Grupo de ligamiento #3 aparece con 5 homólogos, lo que es un error
# Usar función 'overviewSNlinks' para revisar LG=3 en más detalle
# para distinguir homólogos de fragmentos
overviewSNlinks(linkage_df = SN_SN_P2,
                LG_hom_stack = LGHomDf_P2,
                LG = 3,
                LOD_threshold = 0)

# La gráfica anterior mostró que fragmentos #4 y #5 pertenecen al mismo homólogo
# Procediendo a unirlos con la función 'merge_homologues'
LGHomDf_P2_1 <- merge_homologues(LG_hom_stack = LGHomDf_P2,
                                 ploidy = 4,
                                 LG = 3,
                                 mergeList = list(c(4,5)))
                                 # lista de fragmentos a unir: 4 y 5

#Añadiendo marcadores simplex x simplex (tipo 1x1) a los grupos de ligamiento
#ya identificados en...

#Padre 1 (seed parent)
SN_SS_P1 <- linkage(dosage_matrix = filtered_data, 
                    markertype1 = c(1,0),
                    markertype2 = c(1,1),
                    which_parent = 1,
                    ploidy = 4,
                    pairing = "random")

P1_SxS_Assigned <- assign_linkage_group(linkage_df = SN_SS_P1,
                                        LG_hom_stack = LGHomDf_P1,
                                        SN_colname = "marker_a",
                                        unassigned_marker_name = "marker_b",
                                        phase_considered = "coupling",
                                        LG_number = 5,
                                        LOD_threshold = 3,
                                        ploidy = 4)

head(P1_SxS_Assigned)




#Padre 2 (pollen parent)

SN_SS_P2 <- linkage(dosage_matrix = filtered_data, 
                    markertype1 = c(1,0), # simplex x nulliplex
                    markertype2 = c(1,1), # simplex x simplex
                    which_parent = 2,
                    ploidy = 4,
                    pairing = "random")

P2_SxS_Assigned <- assign_linkage_group(linkage_df = SN_SS_P2,
                                        LG_hom_stack = LGHomDf_P2_1,
                                        SN_colname = "marker_a",
                                        unassigned_marker_name = "marker_b",
                                        phase_considered = "coupling",
                                        LG_number = 5,
                                        LOD_threshold = 3,
                                        ploidy = 4)

# Renombrando grupos de ligamiento en base a marcadores simplex x simplex
# (tipo 1x1). Este tipo de marcadores está presente en ambos P1 y P2
# Esta nueva denominación representa los nombres de concenso
LGHomDf_P2_2 <- consensus_LG_names(modify_LG = LGHomDf_P2_1, 
                                   template_SxS = P1_SxS_Assigned, 
                                   modify_SxS = P2_SxS_Assigned)

# Aplicando nuevos nombres de concenso a los grupos de ligamiento
# ya identificados en el Padre 2
P2_SxS_Assigned_2 <- assign_linkage_group(linkage_df = SN_SS_P2,
                                        LG_hom_stack = LGHomDf_P2_2,
                                        SN_colname = "marker_a",
                                        unassigned_marker_name = "marker_b",
                                        phase_considered = "coupling",
                                        LG_number = 5,
                                        LOD_threshold = 3,
                                        ploidy = 4)

# Añadiendo marcadores duplex x nulliplex (tipo 2x0) a los grupos de
# ligamiento de concenso

# Padre 1
P1_DxN_Assigned <- assign_linkage_group(linkage_df = SN_DN_P1,
                                        LG_hom_stack = LGHomDf_P1,
                                        SN_colname = "marker_a",
                                        unassigned_marker_name = "marker_b",
                                        phase_considered = "coupling",
                                        LG_number = 5,
                                        LOD_threshold = 3,
                                        ploidy = 4)

# Padre 2
P2_DxN_Assigned <- assign_linkage_group(linkage_df = SN_DN_P2,
                                        LG_hom_stack = LGHomDf_P2_2,
                                        SN_colname = "marker_a",
                                        unassigned_marker_name = "marker_b",
                                        phase_considered = "coupling",
                                        LG_number = 5,
                                        LOD_threshold = 3,
                                        ploidy = 4)

# Añadiendo todos los otros marcadores

# Padre 1
marker_assignments_P1 <- homologue_lg_assignment(dosage_matrix = filtered_data,
                                                 assigned_list = list(P1_SxS_Assigned, 
                                                                      P1_DxN_Assigned),
                                                 assigned_markertypes = list(c(1,1), c(2,0)),
                                                 LG_hom_stack = LGHomDf_P1,
                                                 which_parent = 1,
                                                 ploidy = 4,
                                                 pairing = "random",
                                                 convert_palindrome_markers = FALSE,
                                                 LG_number = 5,
                                                 LOD_threshold = 3,
                                                 write_intermediate_files = FALSE)

# data(marker_assignments_P1) 
# head(marker_assignments_P1)




# Padre 2
data(marker_assignments_P2)
head(marker_assignments_P2)



# # Control de calidad: chequear consistencia de marcadores entre padres
# # Marcadores presentes en ambos padres deben ser asignados al mismo grupo de
# # ligamiento en ambos casos. Si no, se descartan
marker_assignments <- check_marker_assignment(marker_assignments_P1,marker_assignments_P2)

class(marker_assignments)

head(marker_assignments)


# Cálculo de ligamiento entre todas las combinaciones de tipos de marcadores
# dentro de cada LG

# Padre 1
all_linkages_list_P1 <- finish_linkage_analysis(marker_assignment = marker_assignments$P1,
                                                dosage_matrix = filtered_data,
                                                which_parent = 1,
                                                convert_palindrome_markers = FALSE,
                                                ploidy = 4,
                                                pairing = "random",
                                                LG_number = 5) 

# Padre 2
# all_linkages_list_P2 <- finish_linkage_analysis(marker_assignment = marker_assignments$P2,
#                                                 dosage_matrix = filtered_data,
#                                                 which_parent = 2,
#                                                 convert_palindrome_markers = TRUE,
#                                                 ploidy = 4,
#                                                 pairing = "random",
#                                                 LG_number = 5)

#Lista de los marcadores en 5 grupos de ligamiento del Padre 1

class(all_linkages_list_P1)

# Vista resumida a la estructura interna del objeto 'all_linkages_list_P1'
# (Lista de 5 data drames)
str(all_linkages_list_P1)

# Creando un mapa integrado de todos los LG
# linkages <- list()
# for(lg in names(all_linkages_list_P1)){
#   linkages[[lg]] <- rbind(all_linkages_list_P1[[lg]], all_linkages_list_P2[[lg]])
# }
# 
# integrated.maplist <- MDSMap_from_list(linkages)

data(integrated.maplist)
summary(integrated.maplist)

#Graficando los 5 grupos de linkage
plot_map(maplist = integrated.maplist)











