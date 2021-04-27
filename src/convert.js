
import simplify from './simplify.js';
import createFeature from './feature.js';

// converts GeoJSON feature into an intermediate projected JSON vector format with simplification data

export default function convert(data, options) {
    const features = [];
    if (data.type === 'FeatureCollection') {
        for (let i = 0; i < data.features.length; i++) {
            convertFeature(features, data.features[i], options, i);
        }

    } else if (data.type === 'Feature') {
        convertFeature(features, data, options);

    } else {
        // single geometry or a geometry collection
        convertFeature(features, {geometry: data}, options);
    }

    return features;
}

/**
 * 转换要素，主要是坐标投影变换
 * @param {*} features 传回外部的数据
 * @param {*} geojson 要处理的geojson数据
 * @param {*} options 参数
 * @param {*} options.crs 投影参数， 要区分4326以及3857
 * @param {*} index 处理的geojson数据的序号
 */
function convertFeature(features, geojson, options, index) {
    if (!geojson.geometry) return;

    const crs = options.crs;
    const coords = geojson.geometry.coordinates;
    const type = geojson.geometry.type;
    const tolerance = Math.pow(options.tolerance / ((1 << options.maxZoom) * options.extent), 2);
    let geometry = [];
    let id = geojson.id;
    if (options.promoteId) {
        id = geojson.properties[options.promoteId];
    } else if (options.generateId) {
        id = index || 0;
    }
    if (type === 'Point') {
        convertPoint(coords, geometry, crs);

    } else if (type === 'MultiPoint') {
        for (const p of coords) {
            convertPoint(p, geometry, crs);
        }

    } else if (type === 'LineString') {
        convertLine(coords, geometry, tolerance, false, crs);

    } else if (type === 'MultiLineString') {
        if (options.lineMetrics) {
            // explode into linestrings to be able to track metrics
            for (const line of coords) {
                geometry = [];
                convertLine(line, geometry, tolerance, false, crs);
                features.push(createFeature(id, 'LineString', geometry, geojson.properties));
            }
            return;
        } else {
            convertLines(coords, geometry, tolerance, false, crs);
        }

    } else if (type === 'Polygon') {
        convertLines(coords, geometry, tolerance, true, crs);

    } else if (type === 'MultiPolygon') {
        for (const polygon of coords) {
            const newPolygon = [];
            convertLines(polygon, newPolygon, tolerance, true, crs);
            geometry.push(newPolygon);
        }
    } else if (type === 'GeometryCollection') {
        for (const singleGeometry of geojson.geometry.geometries) {
            convertFeature(features, {
                id,
                geometry: singleGeometry,
                properties: geojson.properties
            }, options, index);
        }
        return;
    } else {
        throw new Error('Input data is not a valid GeoJSON object.');
    }

    features.push(createFeature(id, type, geometry, geojson.properties));
}

function convertPoint(coords, out, crs) {
    out.push(projectX(coords[0], crs), projectY(coords[1], crs), 0);
}

function convertLine(ring, out, tolerance, isPolygon, crs) {
    let x0, y0;
    let size = 0;

    for (let j = 0; j < ring.length; j++) {
        const x = projectX(ring[j][0], crs);
        const y = projectY(ring[j][1], crs);

        out.push(x, y, 0);

        if (j > 0) {
            if (isPolygon) {
                size += (x0 * y - x * y0) / 2; // area
            } else {
                size += Math.sqrt(Math.pow(x - x0, 2) + Math.pow(y - y0, 2)); // length
            }
        }
        x0 = x;
        y0 = y;
    }

    const last = out.length - 3;
    out[2] = 1;
    simplify(out, 0, last, tolerance);
    out[last + 2] = 1;

    out.size = Math.abs(size);
    out.start = 0;
    out.end = out.size;
}

function convertLines(rings, out, tolerance, isPolygon, crs) {
    for (let i = 0; i < rings.length; i++) {
        const geom = [];
        convertLine(rings[i], geom, tolerance, isPolygon, crs);
        out.push(geom);
    }
}

/**
 * 针对X方向投影
 * @param {*} x 
 * @param {*} crs
 * @returns 返回投影后的坐标，注意这里的坐标不是平面坐标值，而是矩阵值，即 [0 ~ 1]之间
 */
function projectX(x, crs) {
    let transformX;
    if (!crs || crs === 'EPSG:3857'){
        transformX = x / 360 + 0.5;
    } else {
        transformX = x / 360 + 0.5;
    }
    return transformX;
}

/**
 * 针对y方向投影
 * @param {*} y
 * @param {*} crs
 * @returns 返回投影后的坐标，注意这里的坐标不是平面坐标值，而是矩阵值，即 [0 ~ 1]之间
 */
function projectY(y, crs) {
    let transformY;
    if (!crs || crs === 'EPSG:3857'){
        const sin = Math.sin(y * Math.PI / 180);
        transformY = 0.5 - 0.25 * Math.log((1 + sin) / (1 - sin)) / Math.PI;    
    } else {
        transformY = (90 - y) / 360;
    }
    return transformY < 0 ? 0 : transformY > 1 ? 1 : transformY;
}
