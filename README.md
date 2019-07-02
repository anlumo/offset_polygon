# offset-polygon

An implementation of the algorithm explained in [Xiaorui Chen, Sara McMains: Polygon Offsetting Using Winding Numbers](https://mcmains.me.berkeley.edu/pubs/DAC05OffsetPolygon.pdf).

The winding number algorithm was adapted from [this page](http://geomalgorithms.com/a03-_inclusion.html), but read the Notes section.

The code itself was written by Andreas Monitzer <andreas@monitzer.com>.

## What does it do?

It allows you to shrink and expand a polygon, like drawing an outline around it. It is also adding arcs to sharp corners with a parameter to control the number of arc points to add, since it outputs polygons only.

## Dependencies

The crate uses [geo-types](https://crates.io/crates/geo-types) in version 0.4 for its data types. The reason is that the author needs to integrate with [geo-booleanop](https://crates.io/crates/geo-booleanop), but it's not really necessary for the operation itself.

## Notes

There are a few magic numbers in the algorithm right now, including the winding number calculation. Initially I used the value returned by `epsilon()`, but it turns out that this fails for a lot of cases (some of which are included as test cases). I don't know why this is, and this might be a problem for different scales. The values right now are optimized for the scale of pixels on a normal screen.

### How can I help?

Just open up a ticket and/or a pull request on this github project. Make sure you explain what you want to do and why.

## License

Licensed under either of

 * Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)
 * [Mozilla Public License 2.0](https://www.mozilla.org/en-US/MPL/2.0/)

at your option.

Dual MIT/Apache2 is strictly more permissive
