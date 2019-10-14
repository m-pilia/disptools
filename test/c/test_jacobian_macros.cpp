#include "gtest/gtest.h"

#include <vector>

extern "C" {
#include "disptools.h"
#include "jacobian_macros.h"
}

namespace {

class JacobianMacrosTest : public ::testing::Test {
protected:
    template<typename T>
    struct Point {
        T x;
        T y;
        T z;
        T d;
    };

    const Point<size_t> size = {5, 6, 7, 3};
    const Point<FLOATING> spacing = {2.0f, 3.0f, 4.0f, 0.0f};
    std::vector<FLOATING> data;
    Image image;

    // Get a component from a voxel of a vector image
    FLOATING get_component(const Point<size_t> point)
    {
        return data[point.d * size.z * size.y * size.x +
                    point.z * size.y * size.x +
                    point.y * size.x +
                    point.x];
    }

    // Get a second order directional derivative
    FLOATING derivative2(
            const FLOATING spacing,
            const Point<size_t> p,
            const Point<size_t> d1,
            const Point<size_t> d2,
            const size_t component)
    {
        const Point<size_t> p0 = {p.x,        p.y,        p.z,        component};
        const Point<size_t> p1 = {p.x + d1.x, p.y + d1.y, p.z + d1.z, component};
        const Point<size_t> p2 = {p.x + d2.x, p.y + d2.y, p.z + d2.z, component};
        return (1.0f / spacing) *
            (-1.5f * get_component(p0) +
             +2.0f * get_component(p1) +
             -0.5f * get_component(p2));
    }

    void SetUp() override
    {
        const size_t data_size = size.z * size.y * size.x * size.d;
        data.reserve(data_size);
        for (int i = 0; i < data_size; ++i) {
            data[i] = static_cast<FLOATING>(i * i);
        }

        image = {
            size.d,
            size.x, size.y, size.z,
            spacing.x, spacing.y, spacing.z,
            data.data(),
        };
    }
};

TEST_F(JacobianMacrosTest, dfx_dx_2f)
{
    const Point<size_t> p = {2, 1, 3, 0};
    FLOATING result = dfx_dx_2f(image, p.x, p.y, p.z, 1.0f / spacing.x);
    FLOATING expected = derivative2(spacing.x, p, {1, 0, 0}, {2, 0, 0}, 0);

    EXPECT_NEAR(result, expected, 1e-6);
}

TEST_F(JacobianMacrosTest, dfx_dy_2f)
{
    const Point<size_t> p = {2, 1, 3, 0};
    FLOATING result = dfx_dy_2f(image, p.x, p.y, p.z, 1.0f / spacing.x);
    FLOATING expected = derivative2(spacing.x, p, {0, 1, 0}, {0, 2, 0}, 0);

    EXPECT_NEAR(result, expected, 1e-6);
}

TEST_F(JacobianMacrosTest, dfx_dz_2f)
{
    const Point<size_t> p = {2, 1, 3, 0};
    FLOATING result = dfx_dz_2f(image, p.x, p.y, p.z, 1.0f / spacing.x);
    FLOATING expected = derivative2(spacing.x, p, {0, 0, 1}, {0, 0, 2}, 0);

    EXPECT_NEAR(result, expected, 1e-6);
}

TEST_F(JacobianMacrosTest, dfy_dx_2f)
{
    const Point<size_t> p = {2, 1, 3, 0};
    FLOATING result = dfy_dx_2f(image, p.x, p.y, p.z, 1.0f / spacing.y);
    FLOATING expected = derivative2(spacing.y, p, {1, 0, 0}, {2, 0, 0}, 1);

    EXPECT_NEAR(result, expected, 1e-6);
}

TEST_F(JacobianMacrosTest, dfy_dy_2f)
{
    const Point<size_t> p = {2, 1, 3, 0};
    FLOATING result = dfy_dy_2f(image, p.x, p.y, p.z, 1.0f / spacing.y);
    FLOATING expected = derivative2(spacing.y, p, {0, 1, 0}, {0, 2, 0}, 1);

    EXPECT_NEAR(result, expected, 1e-6);
}

TEST_F(JacobianMacrosTest, dfy_dz_2f)
{
    const Point<size_t> p = {2, 1, 3, 0};
    FLOATING result = dfy_dz_2f(image, p.x, p.y, p.z, 1.0f / spacing.y);
    FLOATING expected = derivative2(spacing.y, p, {0, 0, 1}, {0, 0, 2}, 1);

    EXPECT_NEAR(result, expected, 1e-6);
}

TEST_F(JacobianMacrosTest, dfz_dx_2f)
{
    const Point<size_t> p = {2, 1, 3, 0};
    FLOATING result = dfz_dx_2f(image, p.x, p.y, p.z, 1.0f / spacing.z);
    FLOATING expected = derivative2(spacing.z, p, {1, 0, 0}, {2, 0, 0}, 2);

    EXPECT_NEAR(result, expected, 1e-6);
}

TEST_F(JacobianMacrosTest, dfz_dy_2f)
{
    const Point<size_t> p = {2, 1, 3, 0};
    FLOATING result = dfz_dy_2f(image, p.x, p.y, p.z, 1.0f / spacing.z);
    FLOATING expected = derivative2(spacing.z, p, {0, 1, 0}, {0, 2, 0}, 2);

    EXPECT_NEAR(result, expected, 1e-6);
}

TEST_F(JacobianMacrosTest, dfz_dz_2f)
{
    const Point<size_t> p = {2, 1, 3, 0};
    FLOATING result = dfz_dz_2f(image, p.x, p.y, p.z, 1.0f / spacing.z);
    FLOATING expected = derivative2(spacing.z, p, {0, 0, 1}, {0, 0, 2}, 2);

    EXPECT_NEAR(result, expected, 1e-6);
}

} // anonymous namespace
