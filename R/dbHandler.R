#########################################
# Functions to handle Mysql connections #
#########################################

library(RMySQL)
library(DBI)

connect <- function(verbose_=FALSE) {
  # Connect to database
  con <- dbConnect(RMySQL::MySQL(),
                   user="user",
                   password="pwd",
                   dbname="kegg-bottleneck",
                   host="10.7.43.4")

  if (verbose_) {
    if (exists(con)) {
      print(paste0("Connected successfully to db ", dbname_))
    } else {
      print(paste0("Something went wrong while connecting to ", dbname_))
    }
  }

  # return the connection object
  return(con)
}

listTables <- function() {
  # Connect to DB
  con <- connect()

  # Perform the function tasks
  print(dbListTables(con))

  # Close DB connection
  dbDisconnect(con)
}

query <- function(sql_) {
  # Connect to DB
  con <- connect()

  # You can fetch all results:
  res <- dbSendQuery(con, sql_)

  # Store the fetched data
  result <- dbFetch(res)

  # Clear the result
  dbClearResult(res)

  # Close DB connection
  dbDisconnect(con)

  # Return the fetched data
  return(result)
}

insertPathway <- function(data_) {
  # Connect to DB
  con <- connect()

  # Construct the insert statement
  sql <- sprintf("insert into pathway
                 (specie, code, gene_count, timestamp)
                 values ('%s', '%s', %s, NOW());",
                 data_["specie"], data_["code"], data_["gene_count"])

  # Debug purpose
  # print(sql)

  # Send the query
  rs <- dbSendQuery(con, sql)

  # Clear the result
  dbClearResult(rs)

  # Get the last inserted ID
  id <- dbGetQuery(con, "select last_insert_id();")[1,1]

  # Close DB connection
  dbDisconnect(con)

  # Return the last inserted ID
  return(id)
}

insertGene <- function(data_) {
  # Connect to DB
  con <- connect()

  # Construct the insert statement
  sql <- sprintf("insert into gene
                 (pathway_id, name, ko, entrez, description, is_bottleneck, belong_to_specie)
                 values (%s, '%s', '%s', '%s', '%s', %s, %s);",
                 data_["pathway_id"], data_["name"], data_["ko"], data_["entrez"],
                 data_["description"], data_["is_bottleneck"], data_["belong_to_specie"])

  # Debug purpose
  print(sql)

  # Send the query
  rs <- dbSendQuery(con, sql)

  # Clear the result
  dbClearResult(rs)

  # Get the last inserted ID
  id <- dbGetQuery(con, "select last_insert_id();")[1,1]

  # Close DB connection
  dbDisconnect(con)

  # Return the last inserted ID
  return(id)
}
