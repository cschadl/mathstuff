#include <array>
#include <functional>
#include <utility>
#include <type_traits>

namespace detail_
{
	// maybe move this to stlutil, since I'm using it in a couple of places...
	template <typename U, size_t... Is>
	constexpr std::array<U, sizeof...(Is)> make_array_val(U const& value, std::index_sequence<Is...>)
	{
		// The static_cast<void> is in case U has overloaded operator, (who does this?)
		// Also it prevents an annoying compiler warning in GCC
		// warning: left operand of comma operator has no effect [-Wunused-value]
		return {{(static_cast<void>(Is), value)...}};
	}

	// Ganked from https://en.cppreference.com/w/cpp/experimental/make_array
	template<class> struct is_ref_wrapper : std::false_type {};
	template<class T> struct is_ref_wrapper<std::reference_wrapper<T>> : std::true_type {};

	template<class B>
	struct negation : std::integral_constant<bool, !bool(B::value)> { };

	template<class T>
	using not_ref_wrapper = negation<is_ref_wrapper<std::decay_t<T>>>;

	template <class D, class...> struct return_type_helper { using type = D; };

	template<class...> struct conjunction : std::true_type { };
	template<class B1> struct conjunction<B1> : B1 { };
	template<class B1, class... Bn>
	struct conjunction<B1, Bn...>
	    : std::conditional_t<bool(B1::value), conjunction<Bn...>, B1> {};

	// gcc gets pissy about this for some reason but it's not a big deal
	//template<class... B>
	//inline constexpr bool conjunction_v = conjunction<B...>::value;

	template <class... Types>
	struct return_type_helper<void, Types...> : std::common_type<Types...> {
		static_assert(conjunction<not_ref_wrapper<Types>...>::value, "Types cannot contain reference_wrappers when D is void");
	};

	template <class D, class... Types>
	using return_type = std::array<typename return_type_helper<D, Types...>::type, sizeof...(Types)>;
};

template < class D = void, class... Types>
constexpr detail_::return_type<D, Types...> make_array(Types&&... t) { return {std::forward<Types>(t)... }; }
