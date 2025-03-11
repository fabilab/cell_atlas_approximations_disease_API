const atlasapprox_disease = require("../index");

// test("highest measurement returns valid structure", async () => {
//     const result = await atlasapprox_disease.highest_measurement({
//         feature: "IL6",
//         number: 10
//     });
//     expect(result).toBeDefined();
//     expect(Array.isArray(result)).toBe(true);
//     result.forEach(entry => {
//         expect(entry).toHaveProperty("cell_type");
//         expect(entry).toHaveProperty("cell_count");
//         expect(entry).toHaveProperty("tissue_general");
//         expect(entry).toHaveProperty("disease");
//         expect(entry).toHaveProperty("dataset_id");
//         expect(entry).toHaveProperty("expression");
//         expect(typeof entry["cell_count"]).toBe("number");
//     });
// });

test("metadata returns expected fields and values", async () => {
    const result = await atlasapprox_disease.metadata({});

    expect(result).toBeDefined();
    expect(Array.isArray(result)).toBe(true);
    result.forEach(entry => {
        expect(entry).toHaveProperty("unique_id");
        expect(entry).toHaveProperty("dataset_id");
        expect(entry).toHaveProperty("cell_type");
        expect(entry).toHaveProperty("tissue_general");
        expect(entry).toHaveProperty("disease");
        expect(entry).toHaveProperty("development_stage_general");
        expect(entry).toHaveProperty("sex");
        expect(entry).toHaveProperty("cell_count");
    });
});

test("metadata with filters", async () => {
    const result = await atlasapprox_disease.metadata({
        disease: "cancer",
        cell_type: "T cell",
        sex: "male"
    });
    expect(result).toBeDefined();
    expect(Array.isArray(result)).toBe(true);
    result.forEach(entry => {
        expect(entry.disease.toLowerCase()).toContain("cancer"); 
        expect(entry.cell_type.toLowerCase()).toContain("t cell");
        expect(entry.sex.toLowerCase()).toBe("male");
    });
});


