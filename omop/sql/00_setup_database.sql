-- One-time database setup. Run as the postgres superuser:
--   psql -U postgres -f omop/sql/00_setup_database.sql
-- The tablespace keeps the database on M: (the default data folder is on C:).
-- M:/PostgresData must exist, be empty, and be writable by the Postgres service account.

CREATE TABLESPACE m_space LOCATION 'M:/PostgresData';
CREATE DATABASE synthea_omop TABLESPACE m_space;

\connect synthea_omop

CREATE SCHEMA native;   -- raw Synthea tables
CREATE SCHEMA cdm;      -- OMOP CDM + vocabulary
CREATE SCHEMA results;  -- Achilles, DQD, cohorts
CREATE SCHEMA scratch;  -- temporary tables
