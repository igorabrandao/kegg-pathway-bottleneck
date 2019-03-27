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

  # Disconnect form DB
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

  # Disconnect from the database
  dbDisconnect(con)

  # Return the fetched data
  return(result)
}
