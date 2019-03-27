-- Truncate table
SET FOREIGN_KEY_CHECKS = 0; 
TRUNCATE table `kegg-bottleneck`.pathway; 
SET FOREIGN_KEY_CHECKS = 1;

DELETE FROM `kegg-bottleneck`.pathway WHERE `id` >= 89;

-- Select all table data
SELECT * FROM `kegg-bottleneck`.pathway;