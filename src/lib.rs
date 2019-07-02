use geo_types::{LineString, Coordinate};
use num_traits::{Num, NumCast, float::{Float, FloatConst}, FromPrimitive};
use std::ops::{AddAssign, SubAssign};

mod error;
mod intersect;
use intersect::intersect;

#[derive(Debug, Clone, Copy)]
struct Normal<N: Num> {
    x: N,
    y: N,
}

#[derive(Debug, Clone)]
struct Segment<N: Num + Copy + NumCast + PartialOrd> {
    p0: Coordinate<N>,
    p1: Coordinate<N>,
    p1_orig: Coordinate<N>,
    normal: Normal<N>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Index {
    Intersection(usize),
    Connected(usize),
}

fn is_left<N>(p0: Coordinate<N>, p1: Coordinate<N>, p2: Coordinate<N>) -> N
        where N: Num + Copy + NumCast + PartialOrd {
    (p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y)
}

// based on http://geomalgorithms.com/a03-_inclusion.html
fn winding_number<N>(pt: Coordinate<N>, polygon: &Vec<Coordinate<N>>) -> isize
        where N: Num + Copy + NumCast + PartialOrd + Float + FloatConst + FromPrimitive + AddAssign + SubAssign + std::fmt::Debug {
    let mut wn = 0;
    let epsilon = N::from_f32(-0.00001).unwrap(); // oh my
    for idx in 0..polygon.len()-1 {
        let p0 = polygon[idx];
        let p1 = polygon[idx+1];
        if p0.y <= pt.y { // start y <= point.y
            if p1.y > pt.y { // an upwards crossing
                if is_left(p0, p1, pt) >= epsilon { // point left of edge
                    wn += 1; // have a valid up intersect
                }
            }
        } else { // start y > point.y
            if p1.y <= pt.y { // a downward crossing
                if is_left(p0, p1, pt) < epsilon { // point right of edge
                    wn -= 1;
                }
            }
        }
    }
    wn
}

pub fn offset_polygon<N>(polygon: &LineString<N>, offset: N, arcdetail: N) -> Result<Vec<LineString<N>>, error::CombinatorialExplosionError>
        where N: Num + Copy + NumCast + PartialOrd + Float + FloatConst + FromPrimitive + AddAssign + SubAssign + std::fmt::Debug {
    if polygon.0.len() == 0 {
        return Ok(vec![LineString(Vec::new())]);
    }
    let arcstep = N::from_f32(2.0).unwrap() * N::PI() / arcdetail;

    let lines: Vec<Segment<N>> = (0..(polygon.0.len()-1)).map(|idx| {
        let (p0, p1) = (polygon.0[idx], polygon.0[idx+1]);
        let len = ((p0.x - p1.x)*(p0.x - p1.x) + (p0.y - p1.y)*(p0.y - p1.y)).sqrt();
        if len < N::epsilon() {
            None
        } else {
            let normal = Normal {
                x: (p1.y - p0.y) / len,
                y: (p0.x - p1.x) / len,
            };
            Some(Segment {
                p0: Coordinate {
                    x: p0.x + offset * normal.x,
                    y: p0.y + offset * normal.y,
                },
                p1: Coordinate {
                    x: p1.x + offset * normal.x,
                    y: p1.y + offset * normal.y,
                },
                p1_orig: p1,
                normal,
            })
        }
    }).flatten().collect();

    let mut connected = Vec::new();

    for idx in 0..lines.len() {
        let (line0, line1) = (&lines[idx], &lines[(idx+1) % lines.len()]);
        connected.extend_from_slice(&[line0.p0, line0.p1]);
        let startangle = line0.normal.y.atan2(line0.normal.x);
        let mut endangle = line1.normal.y.atan2(line1.normal.x);
        let mut angle = startangle - endangle;
        if angle.is_sign_negative() {
            angle += N::from_f32(2.0).unwrap() * N::PI();
        }
        if offset.is_sign_negative() {
            angle = N::from_f32(2.0).unwrap() * N::PI() - angle;
        }
        if angle < N::PI() { // normals facing outwards
            if angle > N::epsilon() {
                connected.push(line0.p1_orig);
            } else { // lines are facing the same direction, remove one of the two coincident points
                connected.pop();
            }
        } else if angle > N::PI() { // normals facing inwards, add arc
            if offset.is_sign_negative() {
                if endangle > startangle {
                    endangle -= N::from_f32(2.0).unwrap() * N::PI();
                }
                for step in 1..<usize as NumCast>::from(((startangle - endangle)/arcstep).ceil()).unwrap() {
                    let angle = startangle - N::from(step).unwrap() * arcstep;
                    connected.push(Coordinate {
                        x: line0.p1_orig.x + offset * angle.cos(),
                        y: line0.p1_orig.y + offset * angle.sin(),
                    });
                }
            } else {
                if endangle < startangle {
                    endangle += N::from_f32(2.0).unwrap() * N::PI();
                }
                for step in 1..<usize as NumCast>::from(((endangle - startangle)/arcstep).ceil()).unwrap() {
                    let angle = startangle + N::from(step).unwrap() * arcstep;
                    connected.push(Coordinate {
                        x: line0.p1_orig.x + offset * angle.cos(),
                        y: line0.p1_orig.y + offset * angle.sin(),
                    });
                }
            }
        }
    }

    connected.push(connected[0]);

    // find intersections and add them to the indices array
    let mut intersections: Vec<Coordinate<N>> = Vec::new();
    let lookup = |idx: Index, intersections: &Vec<Coordinate<N>>, connected: &Vec<Coordinate<N>>| {
        match idx {
            Index::Intersection(idx) => intersections[idx],
            Index::Connected(idx) => connected[idx],
        }
    };
    let mut indices: Vec<Index> = (0..(connected.len()-1)).map(|idx| Index::Connected(idx)).collect();
    let mut indices_idx = 0;
    while indices_idx < indices.len()-1 {
        let mut p0 = lookup(indices[indices_idx], &intersections, &connected);
        let p1 = lookup(indices[(indices_idx+1) % indices.len()], &intersections, &connected);

        loop {
            // exclude line itself from intersection test
            let mut rest = Vec::new();
            if indices_idx+2 < indices.len() {
                rest.extend_from_slice(&indices[(indices_idx+2)..]);
                rest.extend_from_slice(&indices[0..indices_idx]);
            } else {
                rest.extend_from_slice(&indices[((indices_idx+2) % indices.len())..indices_idx]);
            }
            let rest_points = rest.into_iter().map(|idx| lookup(idx, &intersections, &connected)).collect();
            if let Some(int) = intersect(p0, p1, &rest_points, true) {
                intersections.push(int.point);
                if intersections.len() > 3000 {
                    return Err(error::CombinatorialExplosionError);
                }
                indices.insert(indices_idx+1, Index::Intersection(intersections.len()-1));
                let mut other_indices_idx = (indices_idx+3+int.index) % indices.len();
                if other_indices_idx > indices_idx {
                    other_indices_idx += 1;
                }
                indices.insert(other_indices_idx, Index::Intersection(intersections.len()-1));
                p0 = int.point;
                indices_idx += 1;
            } else {
                break;
            }
        }
        indices_idx += 1;
    }

    // find all regions in this polygon
    let mut regions = Vec::new();
    let mut remaining: Vec<usize> = (0..indices.len()).collect();
    while remaining.len() > 0 {
        let mut indices_idx = remaining[0];

        let mut current_region = Vec::new();
        let start_idx = indices_idx;
        loop {
            let idx = indices[indices_idx];
            if let Some(remaining_idx) = remaining.iter().position(|idx| *idx == indices_idx) {
                remaining.remove(remaining_idx);
            }
            match idx {
                Index::Intersection(idx_num) => {
                    current_region.push(intersections[idx_num]);
                    indices_idx = indices.iter().skip(indices_idx+1).position(|i| *i == idx).map(|i| i + indices_idx + 1).unwrap_or_else(|| {
                         // matching entry is before our entry in the array
                        indices.iter().position(|i| *i == idx).unwrap()
                    });
                },
                Index::Connected(idx_num) => {
                    current_region.push(connected[idx_num]);
                },
            }
            indices_idx = (indices_idx+1) % indices.len();
            if start_idx == indices_idx {
                break;
            }
        }
        if current_region.len() > 0 {
            current_region.push(current_region[0]); // line string has to be closed
            regions.push(current_region);
        }
    }

    connected = indices.into_iter().map(|idx| lookup(idx, &intersections, &connected)).collect();
    connected.push(connected[0]); // line string has to be closed

    let epsilon = N::from_f32(0.01).unwrap();
    Ok(regions.into_iter().filter(|region| {
        for idx in 0..(region.len()-1) {
            let p0 = region[idx];
            let p1 = region[idx+1];
            if (p1.y - p0.y).abs() > epsilon {
                return winding_number(Coordinate { // center point
                    x: (p0.x + p1.x) * N::from_f32(0.5).unwrap(),
                    y: (p0.y + p1.y) * N::from_f32(0.5).unwrap(),
                }, &connected) == 1;
            }
        }
        false
    }).map(|coordinates| LineString(coordinates)).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn it_works() {
        let input = LineString(vec![]);
        assert!(offset_polygon(&input, 0.0, 10.0).is_ok(), "offset_polygon returns not an empty list when there's no input");
    }
    #[test]
    fn triangle_extend() {
        let input = LineString(vec![
            Coordinate {
                x: 0.0,
                y: 0.0,
            },
            Coordinate {
                x: 1.0,
                y: 0.0,
            },
            Coordinate {
                x: 1.0,
                y: 1.0,
            },
            Coordinate {
                x: 0.0,
                y: 0.0,
            },
        ]);
        let result = offset_polygon(&input, 0.1, 10.0).unwrap();
        assert!(result.len() == 1, "Triangle offsetting should result in one polygon");
        let first = result.first().unwrap();
        assert!(first.0.iter().zip(vec![Coordinate { x: 0.0, y: -0.1 }, Coordinate { x: 1.0, y: -0.1 }, Coordinate { x: 1.0587785252292474, y: -0.08090169943749476 }, Coordinate { x: 1.0951056516295155, y: -0.03090169943749474 }, Coordinate { x: 1.1, y: 0.0 }, Coordinate { x: 1.1, y: 1.0 }, Coordinate { x: 1.0809016994374947, y: 1.0587785252292474 }, Coordinate { x: 1.0309016994374947, y: 1.0951056516295155 }, Coordinate { x: 0.9690983005625052, y: 1.0951056516295155 }, Coordinate { x: 0.9292893218813453, y: 1.0707106781186548 }, Coordinate { x: -0.07071067811865475, y: 0.07071067811865475 }, Coordinate { x: -0.09876883405951377, y: 0.0156434465040231 }, Coordinate { x: -0.0891006524188368, y: -0.04539904997395467 }, Coordinate { x: -0.04539904997395469, y: -0.08910065241883679 }]).all(|(p0, p1)| (p0.x - p1.x).abs() < f64::epsilon() && (p0.y - p1.y).abs() < f64::epsilon()), "Incorrect triangle offsetting");
    }
    #[test]
    fn triangle_contract() {
        let input = LineString(vec![
            Coordinate {
                x: 0.0,
                y: 0.0,
            },
            Coordinate {
                x: 1.0,
                y: 0.0,
            },
            Coordinate {
                x: 1.0,
                y: 1.0,
            },
            Coordinate {
                x: 0.0,
                y: 0.0,
            },
        ]);
        let result = offset_polygon(&input, -0.1, 10.0).unwrap();
        assert!(result.len() == 1, "Triangle offsetting should result in one polygon");
        let first = result.first().unwrap();
        assert!(first.0.iter().zip(vec![Coordinate { x: 0.9, y: 0.1 }, Coordinate { x: 0.9, y: 0.7585786437626904 }, Coordinate { x: 0.24142135623730954, y: 0.1 }]).all(|(p0, p1)| (p0.x - p1.x).abs() < f64::epsilon() && (p0.y - p1.y).abs() < f64::epsilon()), "Incorrect triangle offsetting");
    }
    #[test]
    fn rectangle() {
        let input = LineString(vec![
            Coordinate { x: 910.0, y: 840.0 }, Coordinate { x: 1890.0, y: 840.0 }, Coordinate { x: 1890.0, y: 1190.0 }, Coordinate { x: 910.0, y: 1190.0 }, Coordinate { x: 910.0, y: 840.0 },
        ]);
        let result = offset_polygon(&input, 4.0, 20.0).unwrap();
        assert!(result.len() == 1, "rectangle input should result in one polygon");
        let first = result.first().unwrap();
        assert!(first.0.iter().zip(vec![
            Coordinate { x: 910.0, y: 836.0 }, Coordinate { x: 1890.0, y: 836.0 }, Coordinate { x: 1891.2360679774997, y: 836.1957739348194 }, Coordinate { x: 1892.3511410091699, y: 836.7639320225002 }, Coordinate { x: 1893.2360679774997, y: 837.6488589908301 }, Coordinate { x: 1893.8042260651805, y: 838.7639320225002 }, Coordinate { x: 1894.0, y: 840.0 }, Coordinate { x: 1894.0, y: 1190.0 }, Coordinate { x: 1893.8042260651805, y: 1191.2360679774997 }, Coordinate { x: 1893.2360679774997, y: 1192.3511410091699 }, Coordinate { x: 1892.3511410091699, y: 1193.2360679774997 }, Coordinate { x: 1891.2360679774997, y: 1193.8042260651805 }, Coordinate { x: 1890.0, y: 1194.0 }, Coordinate { x: 910.0, y: 1194.0 }, Coordinate { x: 908.7639320225002, y: 1193.8042260651805 }, Coordinate { x: 907.6488589908301, y: 1193.2360679774997 }, Coordinate { x: 906.7639320225002, y: 1192.3511410091699 }, Coordinate { x: 906.1957739348194, y: 1191.2360679774997 }, Coordinate { x: 906.0, y: 1190.0 }, Coordinate { x: 906.0, y: 840.0 }, Coordinate { x: 906.1957739348194, y: 838.7639320225002 }, Coordinate { x: 906.7639320225002, y: 837.6488589908301 }, Coordinate { x: 907.6488589908301, y: 836.7639320225002 }, Coordinate { x: 908.7639320225002, y: 836.1957739348194 }, Coordinate { x: 910.0, y: 836.0 },
        ]).all(|(p0, p1)| (p0.x - p1.x).abs() < f64::epsilon() && (p0.y - p1.y).abs() < f64::epsilon()), "Incorrect offsetting for rectangle");
    }
    #[test]
    fn complex() {
        let input = LineString(vec![
            Coordinate { x: 490.0, y: 210.0 }, Coordinate { x: 1260.0, y: 210.0 }, Coordinate { x: 1260.0, y: 433.2142857142857 }, Coordinate { x: 1260.0, y: 538.2142857142858 }, Coordinate { x: 1260.0, y: 1190.0 }, Coordinate { x: 490.0, y: 1190.0 }, Coordinate { x: 484.2650146484375, y: 1189.7904052734375 }, Coordinate { x: 473.10126876831055, y: 1188.4620895385742 }, Coordinate { x: 462.30468368530273, y: 1186.069725036621 }, Coordinate { x: 451.87491607666016, y: 1182.6432037353516 }, Coordinate { x: 441.8116226196289, y: 1178.212417602539 }, Coordinate { x: 432.1144599914551, y: 1172.807258605957 }, Coordinate { x: 422.78308486938477, y: 1166.457618713379 }, Coordinate { x: 413.81715393066406, y: 1159.1933898925781 }, Coordinate { x: 405.21632385253906, y: 1151.0444641113281 }, Coordinate { x: 396.98025131225586, y: 1142.0407333374023 }, Coordinate { x: 389.10859298706055, y: 1132.2120895385742 }, Coordinate { x: 381.6010055541992, y: 1121.5884246826172 }, Coordinate { x: 374.45714569091797, y: 1110.1996307373047 }, Coordinate { x: 367.6766700744629, y: 1098.0755996704102 }, Coordinate { x: 361.2592353820801, y: 1085.246223449707 }, Coordinate { x: 355.2044982910156, y: 1071.7413940429688 }, Coordinate { x: 349.5121154785156, y: 1057.5910034179688 }, Coordinate { x: 344.1817436218262, y: 1042.8249435424805 }, Coordinate { x: 339.21303939819336, y: 1027.4731063842773 }, Coordinate { x: 334.6056594848633, y: 1011.5653839111328 }, Coordinate { x: 330.35926055908203, y: 995.1316680908203 }, Coordinate { x: 326.4734992980957, y: 978.2018508911133 }, Coordinate { x: 322.9480323791504, y: 960.8058242797852 }, Coordinate { x: 319.7825164794922, y: 942.9734802246094 }, Coordinate { x: 316.9766082763672, y: 924.7347106933594 }, Coordinate { x: 314.5299644470215, y: 906.1194076538086 }, Coordinate { x: 312.4422416687012, y: 887.1574630737305 }, Coordinate { x: 310.71309661865234, y: 867.8787689208984 }, Coordinate { x: 308.7462463378906, y: 838.4624633789062 }, Coordinate { x: 307.435302734375, y: 798.36328125 }, Coordinate { x: 307.552490234375, y: 757.53515625 }, Coordinate { x: 309.0950622558594, y: 716.2172241210938 }, Coordinate { x: 312.0602722167969, y: 674.6486206054688 }, Coordinate { x: 316.44537353515625, y: 633.0684814453125 }, Coordinate { x: 322.24761962890625, y: 591.7159423828125 }, Coordinate { x: 327.5717887878418, y: 561.0149459838867 }, Coordinate { x: 331.53319549560547, y: 540.7335662841797 }, Coordinate { x: 335.8473434448242, y: 520.6436004638672 }, Coordinate { x: 340.51388931274414, y: 500.77494049072266 }, Coordinate { x: 345.5324897766113, y: 481.15747833251953 }, Coordinate { x: 350.9028015136719, y: 461.82110595703125 }, Coordinate { x: 356.6244812011719, y: 442.79571533203125 }, Coordinate { x: 362.6971855163574, y: 424.11119842529297 }, Coordinate { x: 369.1205711364746, y: 405.79744720458984 }, Coordinate { x: 375.89429473876953, y: 387.8843536376953 }, Coordinate { x: 383.0180130004883, y: 370.4018096923828 }, Coordinate { x: 390.49138259887695, y: 353.3797073364258 }, Coordinate { x: 398.31406021118164, y: 336.84793853759766 }, Coordinate { x: 406.48570251464844, y: 320.8363952636719 }, Coordinate { x: 415.00596618652344, y: 305.3749694824219 }, Coordinate { x: 423.87450790405273, y: 290.4935531616211 }, Coordinate { x: 433.0909843444824, y: 276.22203826904297 }, Coordinate { x: 442.6550521850586, y: 262.59031677246094 }, Coordinate { x: 452.56636810302734, y: 249.62828063964844 }, Coordinate { x: 462.82458877563477, y: 237.3658218383789 }, Coordinate { x: 473.42937088012695, y: 225.83283233642578 }, Coordinate { x: 484.38037109375, y: 215.0592041015625 }, Coordinate { x: 490.0, y: 210.0 },
        ]);

        let result = offset_polygon(&input, 4.0, 20.0).unwrap();
        assert!(result.len() == 1, "Complex input should result in one polygon");
        let first = result.first().unwrap();
        assert!(first.0.iter().zip(vec![
            Coordinate { x: 490.0, y: 206.0 }, Coordinate { x: 1260.0, y: 206.0 }, Coordinate { x: 1261.2360679774997, y: 206.19577393481939 }, Coordinate { x: 1262.3511410091699, y: 206.7639320225002 }, Coordinate { x: 1263.2360679774997, y: 207.64885899083012 }, Coordinate { x: 1263.8042260651805, y: 208.7639320225002 }, Coordinate { x: 1264.0, y: 210.0 }, Coordinate { x: 1264.0, y: 433.2142857142857 }, Coordinate { x: 1264.0, y: 538.2142857142858 }, Coordinate { x: 1264.0, y: 1190.0 }, Coordinate { x: 1263.8042260651805, y: 1191.2360679774997 }, Coordinate { x: 1263.2360679774997, y: 1192.3511410091699 }, Coordinate { x: 1262.3511410091699, y: 1193.2360679774997 }, Coordinate { x: 1261.2360679774997, y: 1193.8042260651805 }, Coordinate { x: 1260.0, y: 1194.0 }, Coordinate { x: 490.0, y: 1194.0 }, Coordinate { x: 489.8539107738989, y: 1193.997331352042 }, Coordinate { x: 484.1189254223364, y: 1193.7877366254795 }, Coordinate { x: 483.79240923739854, y: 1193.7623876658276 }, Coordinate { x: 472.6286633572716, y: 1192.4340719309644 }, Coordinate { x: 472.235917415522, y: 1192.3673637973461 }, Coordinate { x: 461.4393323325142, y: 1189.974999295393 }, Coordinate { x: 461.05620284228667, y: 1189.8698955044583 }, Coordinate { x: 450.6264352336441, y: 1186.4433742031888 }, Coordinate { x: 450.263066528529, y: 1186.3040698943623 }, Coordinate { x: 440.19977307149776, y: 1181.8732837615498 }, Coordinate { x: 439.8641392595229, y: 1181.7063135028359 }, Coordinate { x: 430.1669766313491, y: 1176.3011545062539 }, Coordinate { x: 429.8641767978316, y: 1176.1142550602183 }, Coordinate { x: 420.5328016757613, y: 1169.7646151676402 }, Coordinate { x: 420.2650159801387, y: 1169.5655648308352 }, Coordinate { x: 411.299085041418, y: 1162.3013360100344 }, Coordinate { x: 411.06603817092304, y: 1162.0970706733294 }, Coordinate { x: 402.46520809279804, y: 1153.9481448920794 }, Coordinate { x: 402.26487638813836, y: 1153.7442711164458 }, Coordinate { x: 394.02880384785516, y: 1144.74054034252 }, Coordinate { x: 393.8581313386792, y: 1144.541206666713 }, Coordinate { x: 385.9864730134839, y: 1134.7125628678848 }, Coordinate { x: 385.8419576995105, y: 1134.5205725682135 }, Coordinate { x: 378.33437026664916, y: 1123.8969077122565 }, Coordinate { x: 378.2124753965426, y: 1123.7139512928131 }, Coordinate { x: 371.0686155332613, y: 1112.3251573475006 }, Coordinate { x: 370.96601961331993, y: 1112.152075028524 }, Coordinate { x: 364.18554399686485, y: 1100.0280439616295 }, Coordinate { x: 364.09926771076067, y: 1099.865066721847 }, Coordinate { x: 357.68183301837786, y: 1087.0356905011438 }, Coordinate { x: 357.60928364250293, y: 1086.882637913942 }, Coordinate { x: 351.5545465514385, y: 1073.3778085072038 }, Coordinate { x: 351.49351299097424, y: 1073.2342396823021 }, Coordinate { x: 345.80113017847424, y: 1059.0838490573021 }, Coordinate { x: 345.7497520513173, y: 1058.9491718322975 }, Coordinate { x: 340.4193801946278, y: 1044.1831119568092 }, Coordinate { x: 340.37610593921517, y: 1044.0566585373924 }, Coordinate { x: 335.40740171558235, y: 1028.7048213791893 }, Coordinate { x: 335.3709444854647, y: 1028.585898674682 }, Coordinate { x: 330.7635645721346, y: 1012.6781762015376 }, Coordinate { x: 330.7328605582804, y: 1012.5660978933767 }, Coordinate { x: 326.48646163249913, y: 996.1323820730642 }, Coordinate { x: 326.460632920877, y: 996.0264879266665 }, Coordinate { x: 322.5748716598907, y: 979.0966707269595 }, Coordinate { x: 322.5531942790427, y: 978.9963372393836 }, Coordinate { x: 319.0277273600974, y: 961.6003106280555 }, Coordinate { x: 319.0096044039095, y: 961.5049559439516 }, Coordinate { x: 315.8440885042513, y: 943.6726118887758 }, Coordinate { x: 315.8290278535101, y: 943.581697031538 }, Coordinate { x: 313.0231196503851, y: 925.342927500288 }, Coordinate { x: 313.0107156855551, y: 925.2559553174486 }, Coordinate { x: 310.5640718562094, y: 906.6406522778979 }, Coordinate { x: 310.55399054475373, y: 906.5571650392243 }, Coordinate { x: 308.4662677664334, y: 887.5952204591462 }, Coordinate { x: 308.45823455756454, y: 887.5147967251163 }, Coordinate { x: 306.7290895075157, y: 868.2361025722843 }, Coordinate { x: 306.7220079614456, y: 868.145623432413 }, Coordinate { x: 304.7551576806839, y: 838.7293178904208 }, Coordinate { x: 304.74838222858386, y: 838.5931636601482 }, Coordinate { x: 303.43743862506824, y: 798.493981531242 }, Coordinate { x: 303.4353192111048, y: 798.3518002410353 }, Coordinate { x: 303.5525067111048, y: 757.5236752410353 }, Coordinate { x: 303.5552750088276, y: 757.3859234096066 }, Coordinate { x: 305.097847030312, y: 716.0679912807003 }, Coordinate { x: 305.1052003538899, y: 715.9326156124139 }, Coordinate { x: 308.0704103148274, y: 674.3640120967889 }, Coordinate { x: 308.08233259815506, y: 674.2291013681985 }, Coordinate { x: 312.46743391651444, y: 632.6489622080422 }, Coordinate { x: 312.4841762830664, y: 632.5126790253505 }, Coordinate { x: 318.2864223768164, y: 591.1601399628505 }, Coordinate { x: 318.3064450906345, y: 591.0324635909847 }, Coordinate { x: 323.63061424957004, y: 560.3314671920589 }, Coordinate { x: 323.64597436945553, y: 560.2481466930744 }, Coordinate { x: 327.6073810772192, y: 539.9667669933674 }, Coordinate { x: 327.6223513493434, y: 539.8937460274299 }, Coordinate { x: 331.9364992985621, y: 519.8037802071174 }, Coordinate { x: 331.9533066012378, y: 519.7290092651388 }, Coordinate { x: 336.6198524691577, y: 499.8603492919942 }, Coordinate { x: 336.63868702806184, y: 499.78357413422617 }, Coordinate { x: 341.657287491929, y: 480.16611197602305 }, Coordinate { x: 341.6783717367669, y: 480.087069936184 }, Coordinate { x: 347.04868347382745, y: 460.7506975606957 }, Coordinate { x: 347.07227651264924, y: 460.66911714782316 }, Coordinate { x: 352.79395620014924, y: 441.64372652282316 }, Coordinate { x: 352.82035906975545, y: 441.5593275250436 }, Coordinate { x: 358.893063384941, y: 422.87481061830533 }, Coordinate { x: 358.9226251549729, y: 422.7873048954861 }, Coordinate { x: 365.34601077509006, y: 404.473553674783 }, Coordinate { x: 365.3791355731061, y: 404.3826468917677 }, Coordinate { x: 372.15285917540103, y: 386.4695533248732 }, Coordinate { x: 372.1900141745596, y: 386.37494795941774 }, Coordinate { x: 379.31373243627837, y: 368.89240401410524 }, Coordinate { x: 379.35545679752823, y: 368.7938040766715 }, Coordinate { x: 386.8288263959169, y: 351.7717017207145 }, Coordinate { x: 386.87574107038193, y: 351.6688197488204 }, Coordinate { x: 394.6984186826866, y: 335.1370509499923 }, Coordinate { x: 394.75123683983315, y: 335.02961798704195 }, Coordinate { x: 402.92287914329995, y: 319.01807471311616 }, Coordinate { x: 402.9824167731871, y: 318.90585414427727 }, Coordinate { x: 411.5026804450621, y: 303.44442836302727 }, Coordinate { x: 411.56986635888444, y: 303.3272346425319 }, Coordinate { x: 420.43840807641374, y: 288.4458183217311 }, Coordinate { x: 420.5142920642168, y: 288.32354175388247 }, Coordinate { x: 429.7307685046465, y: 274.05202686130434 }, Coordinate { x: 429.8165257684177, y: 273.9246656319559 }, Coordinate { x: 439.3805936089939, y: 260.2929441353739 }, Coordinate { x: 439.4775227518404, y: 260.1606446338726 }, Coordinate { x: 449.38883866980916, y: 247.19860850106014 }, Coordinate { x: 449.4983495026084, y: 247.0617144769549 }, Coordinate { x: 459.7565701752158, y: 234.79925567568537 }, Coordinate { x: 459.88015649306186, y: 234.65836573839456 }, Coordinate { x: 470.48493859755405, y: 223.12537623644144 }, Coordinate { x: 470.62413010847314, y: 222.98140739211266 }, Coordinate { x: 481.5751303220962, y: 212.20777915724938 }, Coordinate { x: 481.70406267063237, y: 212.08643212504802 }, Coordinate { x: 487.32369157688237, y: 207.02722802348552 }, Coordinate { x: 488.373316495723, y: 206.34570105534414 }, Coordinate { x: 489.5821725297326, y: 206.02188232890353 }, Coordinate { x: 490.0, y: 206.0 },
        ]).all(|(p0, p1)| (p0.x - p1.x).abs() < f64::epsilon() && (p0.y - p1.y).abs() < f64::epsilon()), "Incorrect offsetting for complex polygon");

        let result = offset_polygon(&input, 40.0, 20.0).unwrap();
        assert!(result.len() == 1, "Complex input should result in one polygon");
    }
    #[test]
    fn complex_contract() {
        let input = LineString(vec![
            Coordinate { x: 490.0, y: 210.0 }, Coordinate { x: 1260.0, y: 210.0 }, Coordinate { x: 1260.0, y: 433.2142857142857 }, Coordinate { x: 1260.0, y: 538.2142857142858 }, Coordinate { x: 1260.0, y: 1190.0 }, Coordinate { x: 490.0, y: 1190.0 }, Coordinate { x: 484.2650146484375, y: 1189.7904052734375 }, Coordinate { x: 473.10126876831055, y: 1188.4620895385742 }, Coordinate { x: 462.30468368530273, y: 1186.069725036621 }, Coordinate { x: 451.87491607666016, y: 1182.6432037353516 }, Coordinate { x: 441.8116226196289, y: 1178.212417602539 }, Coordinate { x: 432.1144599914551, y: 1172.807258605957 }, Coordinate { x: 422.78308486938477, y: 1166.457618713379 }, Coordinate { x: 413.81715393066406, y: 1159.1933898925781 }, Coordinate { x: 405.21632385253906, y: 1151.0444641113281 }, Coordinate { x: 396.98025131225586, y: 1142.0407333374023 }, Coordinate { x: 389.10859298706055, y: 1132.2120895385742 }, Coordinate { x: 381.6010055541992, y: 1121.5884246826172 }, Coordinate { x: 374.45714569091797, y: 1110.1996307373047 }, Coordinate { x: 367.6766700744629, y: 1098.0755996704102 }, Coordinate { x: 361.2592353820801, y: 1085.246223449707 }, Coordinate { x: 355.2044982910156, y: 1071.7413940429688 }, Coordinate { x: 349.5121154785156, y: 1057.5910034179688 }, Coordinate { x: 344.1817436218262, y: 1042.8249435424805 }, Coordinate { x: 339.21303939819336, y: 1027.4731063842773 }, Coordinate { x: 334.6056594848633, y: 1011.5653839111328 }, Coordinate { x: 330.35926055908203, y: 995.1316680908203 }, Coordinate { x: 326.4734992980957, y: 978.2018508911133 }, Coordinate { x: 322.9480323791504, y: 960.8058242797852 }, Coordinate { x: 319.7825164794922, y: 942.9734802246094 }, Coordinate { x: 316.9766082763672, y: 924.7347106933594 }, Coordinate { x: 314.5299644470215, y: 906.1194076538086 }, Coordinate { x: 312.4422416687012, y: 887.1574630737305 }, Coordinate { x: 310.71309661865234, y: 867.8787689208984 }, Coordinate { x: 308.7462463378906, y: 838.4624633789062 }, Coordinate { x: 307.435302734375, y: 798.36328125 }, Coordinate { x: 307.552490234375, y: 757.53515625 }, Coordinate { x: 309.0950622558594, y: 716.2172241210938 }, Coordinate { x: 312.0602722167969, y: 674.6486206054688 }, Coordinate { x: 316.44537353515625, y: 633.0684814453125 }, Coordinate { x: 322.24761962890625, y: 591.7159423828125 }, Coordinate { x: 327.5717887878418, y: 561.0149459838867 }, Coordinate { x: 331.53319549560547, y: 540.7335662841797 }, Coordinate { x: 335.8473434448242, y: 520.6436004638672 }, Coordinate { x: 340.51388931274414, y: 500.77494049072266 }, Coordinate { x: 345.5324897766113, y: 481.15747833251953 }, Coordinate { x: 350.9028015136719, y: 461.82110595703125 }, Coordinate { x: 356.6244812011719, y: 442.79571533203125 }, Coordinate { x: 362.6971855163574, y: 424.11119842529297 }, Coordinate { x: 369.1205711364746, y: 405.79744720458984 }, Coordinate { x: 375.89429473876953, y: 387.8843536376953 }, Coordinate { x: 383.0180130004883, y: 370.4018096923828 }, Coordinate { x: 390.49138259887695, y: 353.3797073364258 }, Coordinate { x: 398.31406021118164, y: 336.84793853759766 }, Coordinate { x: 406.48570251464844, y: 320.8363952636719 }, Coordinate { x: 415.00596618652344, y: 305.3749694824219 }, Coordinate { x: 423.87450790405273, y: 290.4935531616211 }, Coordinate { x: 433.0909843444824, y: 276.22203826904297 }, Coordinate { x: 442.6550521850586, y: 262.59031677246094 }, Coordinate { x: 452.56636810302734, y: 249.62828063964844 }, Coordinate { x: 462.82458877563477, y: 237.3658218383789 }, Coordinate { x: 473.42937088012695, y: 225.83283233642578 }, Coordinate { x: 484.38037109375, y: 215.0592041015625 }, Coordinate { x: 490.0, y: 210.0 },
        ]);
        let result = offset_polygon(&input, -100.0, 20.0).unwrap();
        assert!(result.len() == 1, "Complex input should result in one polygon");
        println!("result: [{}] {:?}", result[0].0.len(), result);
        let first = result.first().unwrap();
        assert!(first.0.iter().zip(vec![
            Coordinate { x: 1160.0, y: 310.0 }, Coordinate { x: 1160.0, y: 433.2142857142857 }, Coordinate { x: 1160.0, y: 433.2142857142857 }, Coordinate { x: 1160.0, y: 538.2142857142858 }, Coordinate { x: 1160.0, y: 538.2142857142858 }, Coordinate { x: 1160.0, y: 1090.0 }, Coordinate { x: 491.95487015907963, y: 1090.0 }, Coordinate { x: 489.8673026705274, y: 1089.7516112627063 }, Coordinate { x: 488.792910658683, y: 1089.5135417906874 }, Coordinate { x: 487.71385124985136, y: 1089.1590353253052 }, Coordinate { x: 486.40119230265503, y: 1088.5810822838698 }, Coordinate { x: 484.69055973565025, y: 1087.627582650634 }, Coordinate { x: 482.4913505001225, y: 1086.1311058172244 }, Coordinate { x: 479.7809439789146, y: 1083.9351249088265 }, Coordinate { x: 476.5894734387631, y: 1080.9113404219213 }, Coordinate { x: 472.9816205961927, y: 1076.9672110438028 }, Coordinate { x: 469.0401531892465, y: 1072.0458491671611 }, Coordinate { x: 464.8534440537063, y: 1066.1214167042547 }, Coordinate { x: 460.50727992896185, y: 1059.1927298855733 }, Coordinate { x: 456.08024714375534, y: 1051.2768422608544 }, Coordinate { x: 451.64167631735677, y: 1042.403499913998 }, Coordinate { x: 447.25121363610583, y: 1032.610762651393 }, Coordinate { x: 442.9593122991198, y: 1021.9417541811705 }, Coordinate { x: 438.8081688067731, y: 1010.4423625800366 }, Coordinate { x: 434.83281253055617, y: 998.1596789161451 }, Coordinate { x: 431.0621852306548, y: 985.1409798366633 }, Coordinate { x: 427.5201291566048, y: 971.4330960608025 }, Coordinate { x: 424.2262511566025, y: 957.0820465050659 }, Coordinate { x: 421.19665724691407, y: 942.1328504166728 }, Coordinate { x: 418.4445657310514, y: 926.6294556027367 }, Coordinate { x: 415.9808128038696, y: 910.6147400212483 }, Coordinate { x: 413.81426614091515, y: 894.1305578135709 }, Coordinate { x: 411.9521613595333, y: 877.2178105688103 }, Coordinate { x: 410.41462503398714, y: 860.0754126142624 }, Coordinate { x: 408.63716155835095, y: 833.4915850598865 }, Coordinate { x: 407.43999339104306, y: 796.87256937449 }, Coordinate { x: 407.5471343036666, y: 759.5446754164375 }, Coordinate { x: 408.96220794896965, y: 721.6417989090759 }, Coordinate { x: 411.6863450617815, y: 683.4527407794883 }, Coordinate { x: 415.7138273744669, y: 645.2635946529438 }, Coordinate { x: 421.05305239344347, y: 607.2110012962893 }, Coordinate { x: 425.920181500002, y: 579.1454555345098 }, Coordinate { x: 429.4999325520822, y: 560.8180542668297 }, Coordinate { x: 433.41730966003814, y: 542.5757563367458 }, Coordinate { x: 437.63886201023803, y: 524.6017364558829 }, Coordinate { x: 442.1590976388982, y: 506.9323579050607 }, Coordinate { x: 446.9716652774589, y: 489.6042018318771 }, Coordinate { x: 452.06923296316756, y: 472.6540714934276 }, Coordinate { x: 457.44334626524164, y: 456.11898174762257 }, Coordinate { x: 463.0842642879559, y: 440.03612833760354 }, Coordinate { x: 468.98077204596916, y: 424.44282967519194 }, Coordinate { x: 475.11996872673336, y: 409.37643150822095 }, Coordinate { x: 481.48703302157094, y: 394.8741619888528 }, Coordinate { x: 488.0649694914805, y: 380.9729212266266 }, Coordinate { x: 494.8343443429819, y: 367.7089854904458 }, Coordinate { x: 501.77302568536385, y: 355.11760208653595 }, Coordinate { x: 508.85595313056723, y: 343.23244718291085 }, Coordinate { x: 516.0549753547518, y: 332.084916611058 }, Coordinate { x: 523.338812698308, y: 321.7032208986425 }, Coordinate { x: 530.673225187708, y: 312.1112635352212 }, Coordinate { x: 532.439413159739, y: 310.0 }, Coordinate { x: 1160.0, y: 310.0 },
        ]).all(|(p0, p1)| (p0.x - p1.x).abs() < f64::epsilon() && (p0.y - p1.y).abs() < f64::epsilon()), "Incorrect offsetting for complex polygon");
    }
}
