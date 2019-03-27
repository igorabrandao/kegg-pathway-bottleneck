-- Truncate table
SET FOREIGN_KEY_CHECKS = 0; 
TRUNCATE table `kegg-bottleneck`.gene; 
SET FOREIGN_KEY_CHECKS = 1;

-- Select all table data
SELECT * FROM `kegg-bottleneck`.gene;

-- Select just genes that belong to the specie pathway
SELECT * FROM `kegg-bottleneck`.gene
WHERE `belong_to_specie` = 1;

-- Select all pathway bottlenecks
SELECT * FROM `kegg-bottleneck`.gene
WHERE `is_bottleneck` = 1;

-- Count how many genes belongs to the specie pathway
SELECT COUNT(`id`) FROM `kegg-bottleneck`.gene
WHERE `belong_to_specie` = 1;

-- Count how many genes are bottlenecks
SELECT COUNT(`id`) FROM `kegg-bottleneck`.gene
WHERE `is_bottleneck` = 1;

DELETE FROM `kegg-bottleneck`.gene WHERE `pathway_id` >= 89;

-- Select all table data
SELECT * FROM `kegg-bottleneck`.gene
WHERE `pathway_id` = 120;