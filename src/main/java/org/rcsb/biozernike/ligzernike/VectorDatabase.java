package org.rcsb.biozernike.ligzernike;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;

public interface VectorDatabase {
    String DB_URL = "jdbc:postgresql://localhost:5432/dalvarez";
    String USER = "dalvarez";
    String PASSWORD = "psqldalvarez";

    // Create the database tables
    void createDatabase();

    void createInfoView();

    // Insert new data
    int insertTarget(String targetName);

    int insertMolecule(String targetName, String moleculeName, boolean active, String smiles);
    int insertMolecule(Integer targetID, String moleculeName, boolean active, String smiles);

    int insertConformer(Integer moleculeID, Integer conformerNumber, double[] geometricVector, String ctab,
            String originalTitle);

    int insertVolumeMap(int conformerID, String mapName);

    // int insertMoments(int mapID, int orderNumber, double[] momentsData);
    int insertMoments(int mapID, int orderNumber, double[] momentsData);

    int insertFPVector(int mapID, double[] fpVectorData);

    // Insert new data
    List<Integer> insertMaps(List<String> mapNames);

    // List<Integer> insertBulkMoments(List<Integer> mapIDs, List<Integer>
    // orderNumbers, List<double[]> momentsData);
    List<Integer> insertBulkMoments(List<Integer> mapIDs, List<Integer> orderNumbers, List<double[]> momentsData);

    // Query data
    void queryData();

    // Close the database connection
    void closeConnection();
}

class VectorDatabaseImpl implements VectorDatabase {

    private Connection connection;

    public VectorDatabaseImpl() {
        try {
            connection = DriverManager.getConnection(DB_URL, USER, PASSWORD);
            createDatabase();
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void createDatabase() {
        try (Statement statement = connection.createStatement()) {
            statement.executeUpdate("CREATE EXTENSION IF NOT EXISTS vector;");

            // Create tables and indexes
            statement.executeUpdate("CREATE TABLE IF NOT EXISTS target (" +
                    "target_id SERIAL PRIMARY KEY," +
                    "target_name VARCHAR(255) NOT NULL UNIQUE)");

            statement.executeUpdate("CREATE TABLE IF NOT EXISTS molecule (" +
                    "molecule_id SERIAL PRIMARY KEY," +
                    "target_id INT NOT NULL REFERENCES target(target_id) ON DELETE CASCADE," +
                    "molecule_name VARCHAR(255) NOT NULL," +
                    "active BOOL NOT NULL," +
                    "smiles TEXT," + // New field for SMILES
                    "UNIQUE (target_id, molecule_name))");

            statement.executeUpdate("CREATE TABLE IF NOT EXISTS conformer (" +
                    "conformer_id SERIAL PRIMARY KEY," +
                    "molecule_id INT NOT NULL REFERENCES molecule(molecule_id) ON DELETE CASCADE," +
                    "conformer_number int NOT NULL," +
                    "original_title VARCHAR(255) NOT NULL," +
                    "geom_vector vector(13)," +
                    "ctab TEXT," + // New field for CTAB
                    "UNIQUE (molecule_id, conformer_number))");

            statement.executeUpdate("CREATE TABLE IF NOT EXISTS volume_map (" +
                    "map_id SERIAL PRIMARY KEY," +
                    "conformer_id INT NOT NULL REFERENCES conformer(conformer_id) ON DELETE CASCADE," +
                    "map_name VARCHAR(255) NOT NULL," +
                    "UNIQUE (conformer_id, map_name))");

            statement.executeUpdate("CREATE TABLE IF NOT EXISTS moments (" +
                    "moment_id SERIAL PRIMARY KEY," +
                    "map_id INT NOT NULL REFERENCES volume_map(map_id) ON DELETE CASCADE," +
                    "order_number INT NOT NULL," +
                    "vector_data vector(444)," +
                    "UNIQUE (map_id, order_number))");

            statement.executeUpdate("CREATE TABLE IF NOT EXISTS fp (" +
                    "fp_vector_id SERIAL PRIMARY KEY," +
                    "map_id INT NOT NULL REFERENCES volume_map(map_id) ON DELETE CASCADE," +
                    "fp_vector vector(72)," +
                    "UNIQUE (map_id))");

            // Add indexes for better performance
            statement.executeUpdate("CREATE INDEX IF NOT EXISTS idx_target_id ON target(target_id)");
            statement.executeUpdate("CREATE INDEX IF NOT EXISTS idx_target_name ON target(target_name)");
            statement.executeUpdate("CREATE INDEX IF NOT EXISTS idx_molecule_id ON molecule(molecule_id)");
            statement.executeUpdate("CREATE INDEX IF NOT EXISTS idx_map_id ON volume_map(map_id)");
            statement.executeUpdate("CREATE INDEX IF NOT EXISTS idx_map_name_conf ON volume_map(map_name, conformer_id)");
            statement.executeUpdate("CREATE INDEX IF NOT EXISTS idx_conformer_id ON conformer(conformer_id)");
            statement.executeUpdate("CREATE INDEX IF NOT EXISTS idx_conformer_num_mol ON conformer(molecule_id, conformer_number)");
            statement.executeUpdate("CREATE INDEX IF NOT EXISTS idx_conformer_orig ON conformer(original_title)");
            statement.executeUpdate(
                    "CREATE INDEX IF NOT EXISTS idx_map_id_order_number ON moments(map_id, order_number)");
            statement.executeUpdate("CREATE INDEX IF NOT EXISTS idx_map_id_fp ON fp(map_id)");

            // Create view
            createInfoView();

        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void createInfoView() {
        try (Statement statement = connection.createStatement()) {
            statement.executeUpdate("CREATE OR REPLACE VIEW info_view AS " +
                    "SELECT " +
                    "    t.target_name, " +
                    "    m.molecule_name, " +
                    "    c.conformer_number, " +
                    "    c.geom_vector, " +
                    "    vm.map_name, " +
                    "    f.fp_vector, " +
                    "    mo.order_number," +
                    "    mo.vector_data AS moments_data " +
                    "FROM " +
                    "    target t " +
                    "JOIN molecule m ON t.target_id = m.target_id " +
                    "JOIN conformer c ON m.molecule_id = c.molecule_id " +
                    "JOIN volume_map vm ON c.conformer_id = vm.conformer_id " +
                    "JOIN moments mo ON vm.map_id = mo.map_id " +
                    "LEFT JOIN fp f ON vm.map_id = f.map_id");
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    @Override
    public int insertTarget(String targetName) {
        int targetId = getTargetId(targetName);

        if (targetId != -1) {
            // Target already exists, return existing target_id
            return targetId;
        }

        try (PreparedStatement preparedStatement = connection.prepareStatement(
                "INSERT INTO target (target_name) VALUES (?)  RETURNING target_id")) {
            preparedStatement.setString(1, targetName.toLowerCase());
            ResultSet resultSet = preparedStatement.executeQuery();
            if (resultSet.next()) {
                targetId = resultSet.getInt(1);
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return targetId;
    }

    public int getTargetId(String targetName) {
        int targetId = -1;
        try (PreparedStatement preparedStatement = connection.prepareStatement(
                "SELECT target_id FROM target WHERE target_name = ?")) {
            preparedStatement.setString(1, targetName.toLowerCase());
            ResultSet resultSet = preparedStatement.executeQuery();
            if (resultSet.next()) {
                targetId = resultSet.getInt(1);
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return targetId;
    }

    public int getMoleculeId(Integer targetID, String moleculeName) {
        int moleculeId = -1;
        try (PreparedStatement preparedStatement = connection.prepareStatement(
                "SELECT molecule_id FROM molecule WHERE molecule_name = ? and target_id = ?;")) {
            preparedStatement.setString(1, moleculeName);
            preparedStatement.setInt(2, targetID);
            ResultSet resultSet = preparedStatement.executeQuery();
            if (resultSet.next()) {
                moleculeId = resultSet.getInt(1);
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return moleculeId;
    }

    public int getConformerId(Integer moleculeID, Integer conformerNumber) {
        int conformerID = -1;
        try (PreparedStatement preparedStatement = connection.prepareStatement(
                "SELECT conformer_id FROM conformer WHERE molecule_id = ? and conformer_number = ?")) {
            preparedStatement.setInt(1, moleculeID);
            preparedStatement.setInt(2, conformerNumber);
            ResultSet resultSet = preparedStatement.executeQuery();
            if (resultSet.next()) {
                conformerID = resultSet.getInt(1);
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return conformerID;
    }

    public int getVolumeId(Integer conformerID, String mapName) {
        int mapID = -1;
        try (PreparedStatement preparedStatement = connection.prepareStatement(
                "SELECT map_id FROM volume_map WHERE conformer_id = ? and map_name = ?")) {
            preparedStatement.setInt(1, conformerID);
            preparedStatement.setString(2, mapName);
            ResultSet resultSet = preparedStatement.executeQuery();
            if (resultSet.next()) {
                mapID = resultSet.getInt(1);
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return mapID;
    }


    @Override
    public int insertMolecule(String targetName, String moleculeName, boolean active, String smiles) {
        int targetID = getTargetId(targetName);
        return insertMolecule(targetID, moleculeName, active, smiles);
    };
    
    @Override
    public int insertMolecule(Integer targetID, String moleculeName, boolean active, String smiles) {
        int moleculeID = getMoleculeId(targetID, moleculeName);

        if (moleculeID != -1) {
            // Target already exists, return existing target_id
            return moleculeID;
        }
        try (PreparedStatement preparedStatement = connection.prepareStatement(
                "INSERT INTO molecule (target_id, molecule_name, active, smiles) VALUES (?, ?, ?, ?) RETURNING molecule_id")) {
            preparedStatement.setInt(1, targetID);
            preparedStatement.setString(2, moleculeName);
            preparedStatement.setBoolean(3, active);
            preparedStatement.setString(4, smiles);
            ResultSet resultSet = preparedStatement.executeQuery();
            if (resultSet.next()) {
                moleculeID = resultSet.getInt(1);
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return moleculeID;
    }

    @Override
    public int insertConformer(Integer moleculeID, Integer conformerNumber, double[] geometricVector, String ctab,
            String originalTitle) {
        int conformerID = getConformerId(moleculeID, conformerNumber);

        if (conformerID != -1) {
            // Target already exists, return existing target_id
            return conformerID;
        }
        try (PreparedStatement preparedStatement = connection.prepareStatement(
                "INSERT INTO conformer (molecule_id, conformer_number, original_title, geom_vector, ctab) VALUES (" +
                        "?, ?, ?, ?, ?)  RETURNING conformer_id")) {
            preparedStatement.setInt(1, moleculeID);
            preparedStatement.setInt(2, conformerNumber);
            preparedStatement.setString(3, originalTitle);
            preparedStatement.setObject(4, geometricVector);
            preparedStatement.setString(5, ctab);
            ResultSet resultSet = preparedStatement.executeQuery();
            if (resultSet.next()) {
                conformerID = resultSet.getInt(1);
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return conformerID;
    }

    @Override
    public int insertVolumeMap(int conformerID, String mapName) {
        int generatedId = -1;
        try (PreparedStatement preparedStatement = connection.prepareStatement(
                "INSERT INTO volume_map (conformer_id, map_name) VALUES (?, ?)  RETURNING map_id")) {
            preparedStatement.setInt(1, conformerID);
            preparedStatement.setString(2, mapName);
            ResultSet resultSet = preparedStatement.executeQuery();
            if (resultSet.next()) {
                generatedId = resultSet.getInt(1);
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return generatedId;
    }


    public int getMomentId(Integer mapID, Integer orderNumber) {
        int generatedId = -1;
        try (PreparedStatement preparedStatement = connection.prepareStatement(
                "SELECT moment_id FROM moments WHERE map_id = ? and order_number = ?")) {
            preparedStatement.setInt(1, mapID);
            preparedStatement.setInt(2, orderNumber);
            ResultSet resultSet = preparedStatement.executeQuery();
            if (resultSet.next()) {
                generatedId = resultSet.getInt(1);
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return generatedId;
    }

    @Override
    public int insertMoments(int mapID, int orderNumber, double[] momentsData) {
        int generatedId = getMomentId(mapID, orderNumber);
        if (generatedId == -1){
            try (PreparedStatement preparedStatement = connection.prepareStatement(
                    "INSERT INTO moments (map_id, order_number, vector_data) VALUES (?, ?, ?) RETURNING moment_id")) {
                preparedStatement.setInt(1, mapID);
                preparedStatement.setInt(2, orderNumber);
                preparedStatement.setObject(3, momentsData);
                //System.out.println("momentsData: "+momentsData[0]+","+momentsData[1]);
                ResultSet resultSet = preparedStatement.executeQuery();
                if (resultSet.next()) {
                    generatedId = resultSet.getInt(1);
                }
            } catch (SQLException e) {
                e.printStackTrace();
            }
        }        
        return generatedId;
    }

    @Override
    public List<Integer> insertMaps(List<String> mapNames) {
        List<Integer> generatedIds = new ArrayList<>();
        try (PreparedStatement preparedStatement = connection.prepareStatement(
                "INSERT INTO volume_map (map_name) VALUES (?) RETURNING map_id")) {
            for (String mapName : mapNames) {
                preparedStatement.setString(1, mapName);
                ResultSet resultSet = preparedStatement.executeQuery();
                if (resultSet.next()) {
                    generatedIds.add(resultSet.getInt(1));
                }
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return generatedIds;
    }

    @Override
    public List<Integer> insertBulkMoments(List<Integer> mapIDs, List<Integer> orderNumbers,
            List<double[]> momentsData) {
        List<Integer> generatedIds = new ArrayList<>();
        try (PreparedStatement preparedStatement = connection.prepareStatement(
                "INSERT INTO moments (map_id, order_number, vector_data) VALUES (?, ?, ?) RETURNING moment_id")) {
            for (int i = 0; i < mapIDs.size(); i++) {
                preparedStatement.setInt(1, mapIDs.get(i));
                preparedStatement.setInt(2, orderNumbers.get(i));
                preparedStatement.setObject(3, momentsData.get(i)); // Use double8 for double
                ResultSet resultSet = preparedStatement.executeQuery();
                if (resultSet.next()) {
                    generatedIds.add(resultSet.getInt(1));
                }
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return generatedIds;
    }

    // private int getMapId(String mapName) {
    // int mapId = -1;
    // try (PreparedStatement preparedStatement = connection.prepareStatement(
    // "SELECT map_id FROM volume_map WHERE map_name = ?")) {
    // preparedStatement.setString(1, mapName);
    // ResultSet resultSet = preparedStatement.executeQuery();
    // if (resultSet.next()) {
    // mapId = resultSet.getInt(1);
    // }
    // } catch (SQLException e) {
    // e.printStackTrace();
    // }
    // return mapId;
    // }

    @Override
    public int insertFPVector(int mapID, double[] fpVectorData) {
        int generatedId = -1;
        try (PreparedStatement preparedStatement = connection.prepareStatement(
                "INSERT INTO fp (map_id, fp_vector) VALUES (?, ?) returning fp_vector_id;")) {
            preparedStatement.setInt(1, mapID);
            preparedStatement.setObject(2, fpVectorData);
            ResultSet resultSet = preparedStatement.executeQuery();
            if (resultSet.next()) {
                generatedId = resultSet.getInt(1);
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return generatedId;
    }

    @Override
    public void queryData() {
        try (Statement statement = connection.createStatement();
                ResultSet resultSet = statement.executeQuery("SELECT * FROM target")) {
            // Process the result set
            while (resultSet.next()) {
                System.out.println("Target ID: " + resultSet.getInt("target_id") +
                        ", Target Name: " + resultSet.getString("target_name"));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void closeConnection() {
        try {
            if (connection != null && !connection.isClosed()) {
                connection.close();
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    public void truncateAll() {
        try (Statement statement = connection.createStatement();
                ResultSet resultSet = statement.executeQuery("TRUNCATE TABLE TARGET CASCADE;")) {
        } catch (SQLException e) {
            e.printStackTrace();
        }

    }

}
